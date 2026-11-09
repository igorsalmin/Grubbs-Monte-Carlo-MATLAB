% Расчет одностороннего критерия Смирнова-Граббса уровня alphaVal
% для размера выборки sampleSizes с построением гистограммы цветной и чб
clc; clear; close all;

alphaVal       = 0.05 * 0.5;  % значение 1-стороннего критерия = 0.5 от 2-стороннего
sampleSizes    = 20;          % одно значение размера выборки n
numSimulations = 1e6;         % количество симуляций Монте-Карло

numSizes = length(sampleSizes);% количество вариантов размеров выборки n 

critValG1 = zeros(numSizes, 1);
critValGn = zeros(numSizes, 1);

for k = 1:numSizes 
    n = sampleSizes(k);
    tstart = tic;
    fprintf('Расчет для n = %d, numSimulations = %.0e,  Время начала: %s\n', n, numSimulations, datetime('now'));

    % генерация numSimulations стандартных нормальных распределений 
    % по n случайных значений в каждой выборке
    % x - это двумерный массив [n, numSimulations]
    x = randn(n, numSimulations);
    
    % среднее значение - одномерный массив [1, numSimulations]
    Mean = mean(x, 1);
    
    % СКО - одномерный массив [1, numSimulations]
    s    = std(x, 0, 1); 

    % Распределение статистик G1 и Gn - одномерные массивы [1, numSimulations]
    localGn = (max(x, [], 1) - Mean) ./ s;
    localG1 = (Mean - min(x, [], 1)) ./ s;
    
    %   Односторонний критерий Смирнова-Граббса уровня alpha для n
    % = двустороннему критерию Смирнова-Граббса уровня alpha*2
    critValG1(k) = quantile(localG1, 1 - alphaVal );
    critValGn(k) = quantile(localGn, 1 - alphaVal );
end

%==========================================================================
% ОФОРМЛЕНИЕ И ЭКСПОРТ РЕЗУЛЬТАТОВ
%==========================================================================

% 1. Создаем красивую таблицу MATLAB
resultsTable = table(sampleSizes', critValG1, critValGn, ...
    'VariableNames', {'SampleSize_n', 'G1_CritVal', 'Gn_CritVal'});

% 2. Выводим результаты в командное окно
fprintf('\n=== РЕЗУЛЬТАТЫ РАСЧЕТА МОНТЕ-КАРЛО (alpha = %.4f) ===\n', alphaVal);
disp(resultsTable);

% 3. Экспортируем в Excel-файл
outputFileName = 'Grubbs_Critical_Values.xlsx';
writetable(resultsTable, outputFileName);
fprintf('Результаты успешно сохранены в файл: %s\n', outputFileName);

%==========================================================================
% ГРАФИКИ
%==========================================================================

% включение поддержки кириллицы
set(0, 'DefaultAxesFontName', 'Times');
set(0, 'DefaultTextFontName', 'Times');

% Монохромная гистограмма
plotHist(localGn',[0.90 0.90 0.90], [0 0 0], [0.2 0.2 0.2]);
% Цветная гистограмма
plotHist(localGn', [102/255 170/255 215/255],  [0.7 0 0.1], 'red');


% =========================================================================
% Гистограмма с красивыми границами и закраской квантиля по форме кривой
% =========================================================================
function plotHist(data, histFaceColor, fitLineColor, quantileFaceColor)
    % 1. Исходные данные
    %data = localGn'; 
    
    % 2. Задаем красивый шаг и рассчитываем границы столбцов
    step = 0.2; 
    xMin = floor(min(data) / step) * step;
    xMax =  ceil(max(data) / step) * step;
    edges = xMin : step : xMax;
    
    % --- РАСЧЕТ КВАНТИЛЯ ---
    q975 = quantile(data, 0.975); 
    
    % 3. Строим правильную гистограмму
    figure;
    hHist = histogram(data, 'BinEdges' , edges, ...
                            'FaceColor', histFaceColor, ... 
                            'FaceAlpha', 1.0, ...
                            'EdgeColor', [0.4 0.4 0.4]);

    % Безопасное изменение толщины линий столбцов для MATLAB версии 2024а
    set(hHist, 'LineWidth', 1.0); 
    hold on;
    
    % 4. Рассчитываем параметры логнормального распределения
    pd = fitdist(data, 'Lognormal');
    
    % --- РАСЧЕТ И ЗАКРАСКА ОБЛАСТИ КВАНТИЛЯ ПО ФОРМЕ КРИВОЙ ---
    % Создаем сетку точек строго ОТ квантиля ДО правого края
    xArea = linspace(q975, xMax, 100);
    % Считаем высоту кривой на этом участке с масштабированием
    yArea = pdf(pd, xArea) * numel(data) * step;
    
    % Закрашиваем область под кривой
    hArea = area(xArea, yArea, 'BaseValue', 0);
    % Настраиваем цвет закраски (тёмно-серый с прозрачностью)
    set(hArea, 'FaceColor', quantileFaceColor, 'FaceAlpha', 0.50, 'EdgeColor', 'none');
    
    % --- НАКЛАДЫВАЕМ ПОЛНУЮ КРИВУЮ РАСПРЕДЕЛЕНИЯ ---
    xGrid = linspace(xMin, xMax, 250);
    yGrid = pdf(pd, xGrid) * numel(data) * step;
    
    % Строим главную линию распределения
    plot(xGrid, yGrid, 'LineWidth', 2.5, 'Color', fitLineColor);
    
    % --- СТРЕЛКА И ПОДПИСЬ КВАНТИЛЯ ---
    yLimits = ylim; 
    
    % --- ВЫНОСКА ДЛЯ MATLAB R2024 ---
    % 1. Определяем координаты начала и конца стрелки (в единицах вашего графика)
    xEnd = q975;
    yEnd = 0;                  % Кончик указывает на ось X в точку квантиля
    
    xStart = q975 + 0.1*0;        % Смещение начала стрелки вправо
    yStart = yLimits(2) * 0.30; % Высота начала стрелки (15% от верха графика)
    
    % Формируем текст подписи
    labelText = sprintf('G_2_0_,  _2_,_5_%% = %.2f', q975);
    
    % 2. Создаем стандартную стрелку
    tg = annotation('textarrow');
    
    % 3. Привязываем её к осям координат графика (gca) и задаем точки
    set(tg, 'Parent', gca, ...
            'X', [xStart, xEnd], ...
            'Y', [yStart, yEnd], ...
            'String', labelText);
    
    % 4. Настройка внешнего вида под ваш стиль
    set(tg, 'FontName', 'Times New Roman', 'FontSize', 14, ...
            'Color', [0.2 0.2 0.2], ...       % Цвет стрелки
            'TextColor', [0 0 0], ...         % Цвет текста
            'LineWidth', 2.0, ...             % Толщина стрелки
            'TextEdgeColor', 'none', ...      % Без рамки вокруг текста
            'TextBackgroundColor', 'none', ...% Прозрачный фон под текстом
            'HorizontalAlignment', 'left');   % Текст пишется вправо от начала стрелки
    
    
    % 5. Оформление осей, сетки и шрифтов
    grid on;
    xlabel('Gn');
    ylabel('Количество элементов');
    
    % Применяем Times New Roman и размер 14 ко всему графику
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
    
    hold off;
    xlim([0, round(xMax)]);
    xticks([0:1:round(xMax)]);
end