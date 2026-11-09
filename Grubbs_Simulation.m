% Расчет одностороннего критерия Смирнова-Граббса уровня alphaVal, который
%      = двустороннему критерию Смирнова-Граббса уровня alphaVal*2
% для различных значений размеров выборок sampleSizes
clear; close all;

alphaVal       = 0.05 * 0.5;           % значение 1-стороннего критерия = 0.5 от 2-стороннего

% sampleSizes - диапазон значений размера выборки n
sampleSizes    = [3:150];
%sampleSizes    = [  3:150,     200:     50:    1000];
%sampleSizes    = [   1500,    2000:   1000:   10000]; 
%sampleSizes    = [  15000,   20000:  10000:  100000]; 
%sampleSizes    = [ 150000,  200000: 100000: 1000000];  
%sampleSizes    = [1500000, 2000000:1000000:10000000]; 


% ДЛЯ ЗАЩИТЫ ОЗУ ВНУТРИ ВОРКЕРОВ:
% Ограничение размера матрицы в ОЗУ одного воркера
maxChunkSize = 50; 

numSubSimulations   = 1e1;  % число случаев, на которые хватит памяти в 1 итерацию
numIteration        = 1e2;  % число итераций, в каждой numSubSimulations 
totalSimulations    = numSubSimulations * numIteration;  % Общее число симуляций

% =========================================================================
% Для экономии ОЗУ будем хранить только «хвосты» распределений G1 и Gn 
% Используем верхнюю границу доверительного интервала биномиального распределения
% =========================================================================
mu_binom  = numSubSimulations * alphaVal;
std_binom = sqrt(numSubSimulations * alphaVal * (1 - alphaVal));

% Берем с запасом в 5.0 сигм (вероятность утери экстремума < 0.00003%)
keepCount = ceil(mu_binom + 3.0 * std_binom); 

% страховка на выход за границы
keepCount = max(1, min(keepCount, numSubSimulations));
% =====================================================================

numSizes  = length(sampleSizes);    % количество вариантов размеров выборки n 
critValG1 = zeros(numSizes, 1);     % выделяем память под критические значения G1
critValGn = zeros(numSizes, 1);     % выделяем память под критические значения Gn
durationTime = zeros(numSizes, 1);

% -------------------------------------------------------------------------
% создаем папку для результатов, если ее еще нет
folderName = 'results';
% Проверяем, существует ли папка, и если нет — создаем её
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
outputFileName = fullfile(folderName, 'Grubbs_Critical_Values.xlsx');
% -------------------------------------------------------------------------

fprintf('Расчет для numSubSimulations = %.0e, numIteration = %.0e\n', numSubSimulations, numIteration);

for k = 1:numSizes 
    n = sampleSizes(k);
    tstart = tic;
    
    if k == 1 
       forcastTimeSecs = n*totalSimulations * 1.8e-9;
    else
       forcastTimeSecs = durationTime(k-1) * sampleSizes(k) / sampleSizes(k-1);
    end
    fprintf('Расчет для n = %7d, Начало: %s.   прогноз финиша %s.      ', n, datetime('now'), datetime('now') + seconds(forcastTimeSecs) );
    
    % Выделяем память только под ТОП самых больших значений
    topG1 = zeros(numIteration, keepCount);
    topGn = zeros(numIteration, keepCount);
    
    % запуск параллельных вычислений
    parfor m = 1:numIteration
        localGn = zeros(numSubSimulations, 1);
        localG1 = zeros(numSubSimulations, 1);
        
        % Внутреннее дробление numSubSimulations на безопасные чанки
        currIdx = 1;
        while currIdx <= numSubSimulations

            % при больших n и totalSimulations генерация numSubSimulations 
            % случайных значений приведет к ошибке памяти, поэтому делим
            % генерации на блоки chunkSize стандартных нормальных распределений 
            % по n случайных значений в каждой выборке
            % результат - двумерный массив [chunkSize, n]
            
            chunkSize = min(maxChunkSize, numSubSimulations - currIdx + 1);
            idxEnd    = currIdx + chunkSize - 1;
                
            % Выделение памяти под маленький безопасный блок
            x = randn(chunkSize, n);
    
            % средние значения - одномерный массив [numSubSimulations, 1]
            mu = mean(x, 2);
    
            % СКО - одномерный массив [numSubSimulations, 1]
            %stdVals  = std(x, 0, 2);
            % Метод смещенных моментов должен быть быстрее чем std
            stdVals = sqrt((sum(x.^2, 2) - n * (mu.^2)) / (n - 1));
    
            localGn(currIdx:idxEnd) = (max(x, [], 2) - mu) ./ stdVals;
            localG1(currIdx:idxEnd) = (mu - min(x, [], 2)) ./ stdVals;
                
            currIdx = idxEnd + 1;
        end

        % Сортируем локальные результаты по возрастанию
        % localSortedGn = sort(localGn);
        % localSortedG1 = sort(localG1);
        % 
        % % Забираем в матрицу ТОЛЬКО хвост (самые большие значения)
        % topGn(m, :) = localSortedGn((end - keepCount + 1):end);
        % topG1(m, :) = localSortedG1((end - keepCount + 1):end);

        % =================================================================
        % ИДЕЯ ОПТИМИЗАЦИИ: ЧАСТИЧНЫЙ ПОИСК ХВОСТОВ ВМЕСТО ПОЛНОЙ СОРТИРОВКИ
        % =================================================================
        % Для Gn: Находим keepCount самых больших элементов.
        % Функция maxk возвращает их в порядке убывания.
        topK_Gn = maxk(localGn, keepCount);
        
        % Для G1: Нам тоже нужны самые большие значения из localG1.
        topK_G1 = maxk(localG1, keepCount);
        
        % Переворачиваем векторы (с помощью flipud), чтобы элементы шли 
        % по возрастанию — в точности как было при старой логике отсечения!
        topGn(m, :) = flipud(topK_Gn);
        topG1(m, :) = flipud(topK_G1);
    end

    % Собираем все топы в один вектор и сортируем еще раз
    topG1 = sort( reshape(topG1, 1, []) ); 
    topGn = sort( reshape(topGn, 1, []) );

    % Находим индекс нужной квантили с конца вектора топ-значений
    targetIdxFromEnd = round(totalSimulations * alphaVal);
    
    critValG1(k) = topG1(end - targetIdxFromEnd + 1);
    critValGn(k) = topGn(end - targetIdxFromEnd + 1);
   
    durationTime(k) = toc(tstart);
    fprintf('Длительность: %7.1f сек.  G1 = %.16f  Gn = %.16f\n', durationTime(k), critValG1(k), critValGn(k));

    % Принудительно очищаем тяжелые переменные в конце шага 'k'
    clear topG1 topGn;

    %==========================================================================
    % 1. Формирование итоговой таблицы результатов
    resultsTable = table(sampleSizes', critValG1, critValGn, durationTime, ...
                        'VariableNames', {'SampleSize_n', 'G1_CritVal', 'Gn_CritVal', 'durationTime'});

    % 3. Экспортируем в Excel-файл
    writetable(resultsTable, outputFileName);
end

%==========================================================================
% 1. Формирование итоговой таблицы результатов
resultsTable = table(sampleSizes', critValG1, critValGn, durationTime, ...
    'VariableNames', {'SampleSize_n', 'G1_CritVal', 'Gn_CritVal', 'durationTime'});

% 2. Выводим результаты в командное окно
fprintf('\n=== КРИТИЧЕСКИХ ЗНАЧЕНИЙ СМИРНОВА-ГРАББСА (ОДНОСТОРОННИЙ) (alpha = %.4f) ===\n', alphaVal);
disp(resultsTable);

% 3. Экспортируем в Excel-файл
writetable(resultsTable, outputFileName);
fprintf('Результаты успешно сохранены в файл: %s\n', outputFileName);

%==========================================================================
% ГРАФИКИ
%==========================================================================

% включение поддержки кириллицы
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

% График значений критерия Смирнова-Граббса в зависимости от n
figure;
subplot(2, 1, 1);
plot(sampleSizes, critValGn);
xlabel("n");
hold on;
grid on;
plot(sampleSizes, critValG1);

subplot(2, 1, 2);
plot(sampleSizes, critValGn - critValG1);
xlabel("n");
grid on;