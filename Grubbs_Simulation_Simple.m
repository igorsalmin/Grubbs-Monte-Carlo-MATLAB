% Расчет одностороннего критерия Смирнова-Граббса уровня alphaVal
% для различных значений размеров выборок sampleSizes
clc; clear; close all;

alphaVal       = 0.05 * 0.5;  % значение 1-стороннего критерия = 0.5 от 2-стороннего
sampleSizes         = 3:150;  % диапазон значений n
numSubSimulations   = 1e4;    % число случаев, на которые хватит памяти в 1 итерацию
numIteration        = 1e2;    % число итераций, в каждой numSubSimulations 
totalSimulations    = numSubSimulations * numIteration;  % Общее число симуляций

% =========================================================================
% Для экономии ОЗУ будем хранить только «хвосты» распределений G1 и Gn 
% Используем верхнюю границу доверительного интервала биномиального распределения
% =========================================================================
mu_binom  = numSubSimulations * alphaVal;
std_binom = sqrt(numSubSimulations * alphaVal * (1 - alphaVal));

% Берем с запасом в 5.0 сигмы (вероятность утери экстремума < 0.00003%)
keepCount = ceil(mu_binom + 5.0 * std_binom); 

% страховка на выход за границы
keepCount = max(1, min(keepCount, numSubSimulations));
% =========================================================================

numSizes  = length(sampleSizes);    % количество вариантов размеров выборки n 
critValG1 = zeros(numSizes, 1);     % выделяем память под критические значения G1
critValGn = zeros(numSizes, 1);     % выделяем память под критические значения Gn
durationTime = zeros(numSizes, 1);

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
    %fprintf('Расчет для n = %d,  Время: %s\n', n, datetime('now'));
    
    % Выделяем память только под ТОП самых больших значений
    topG1 = zeros(numIteration, keepCount);
    topGn = zeros(numIteration, keepCount);
    
    % запуск параллельных вычислений
    parfor m = 1:numIteration
        % генерация numSubSimulations стандартных нормальных распределений 
        % по n случайных значений в каждой выборке
        % результат - двумерный массив [numSubSimulations, n]
        x = randn(numSubSimulations, n);
        % средние значения - одномерный массив [numSubSimulations, 1]
        mu = mean(x, 2);
        % СКО - одномерный массив [numSubSimulations, 1]
        stdVals  = std(x, 0, 2);

        localGn = (max(x, [], 2) - mu  ) ./ stdVals;
        localG1 = (mu   - min(x, [], 2)) ./ stdVals;
       
        % Сортируем локальные результаты по возрастанию
        localSortedGn = sort(localGn);
        localSortedG1 = sort(localG1);
        
        % Забираем в матрицу ТОЛЬКО хвост (самые большие значения)
        topGn(m, :) = localSortedGn((end - keepCount + 1) : end);
        topG1(m, :) = localSortedG1((end - keepCount + 1) : end);
    end
    
    % Собираем все топы в один вектор и сортируем еще раз
    topG1 = topG1(:);
    topGn = topGn(:);
    allTopG1 = sort(topG1);
    allTopGn = sort(topGn);
    
    % Находим индекс нужной квантили с конца вектора топ-значений
    targetIdxFromEnd = round(totalSimulations * alphaVal);
    critValG1(k) = allTopG1(end - targetIdxFromEnd + 1);
    critValGn(k) = allTopGn(end - targetIdxFromEnd + 1);
    
    durationTime(k) = toc(tstart);
    fprintf('Длительность: %7.1f сек.  G1 = %.16f  Gn = %.16f\n', durationTime(k), critValG1(k), critValGn(k));


    % Принудительно очищаем тяжелые переменные в конце шага 'k'
    clear topG1 topGn allTopG1 allTopGn;
end

% 1. Формирование итоговой таблицы результатов
resultsTable = table(sampleSizes', critValG1, critValGn, ...
    'VariableNames', {'SampleSize_n', 'G1_CritVal', 'Gn_CritVal'});

% 2. Выводим результаты в командное окно
fprintf('\n=== КРИТИЧЕСКИЕ ЗНАЧЕНИЯ ОДНОСТОРОННЕГО КРИТЕРИЯ СМИРНОВА-ГРАББСА …  (alpha = %.4f) ===\n', alphaVal);
disp(resultsTable);

% 3. Экспортируем в Excel-файл в подпапку results
folderName = 'results';

% Проверяем, существует ли папка, и если нет — создаем её
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

outputFileName = fullfile(folderName, 'Grubbs_Critical_Values_Simple.xlsx');
writetable(resultsTable, outputFileName);