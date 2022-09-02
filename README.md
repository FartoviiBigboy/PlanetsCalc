# PlanetsCalc
Решение задачи движения N тел в трехмерном пространстве с учетом сил гравитационного притяжения между ними; 
при некоторых упругих соударениях количество тел увеличивается (одно из тел разделяется на 2) с заданной вероятностью, иначе не меняется. 
Массы, радиусы, начальные координаты и векторы скоростей всех тел считать заданными. Вычислить координаты всех тел через время T.
___
## Работа программы:
После запуска программы, на выбор будет 5 пунктов:
1. Choose file from disk – загрузить входные данные с каталога
2. Generate data – сгенерировать входные данные и сохранить их в каталог
3. Calculation from a single- and multi-threaded program – произвести вычисления с загруженными данными (данный пункт запускать ТОЛЬКО после выполнения 1 пункта)
4. Compare results – сравнить полученные результаты для одно- и много-поточной реализации и сохранить полученные результаты вычислений (данный пункт запускать ТОЛЬКО после выполнения 3 пункта)
5. Exit – выйти из программы
___
### 1 пункт:
Предлагается ввести имя файла, с которого будут считываться входные данные, файл должен находиться в том же каталоге, что и исполняемая программа.
___
### 2 пункт:
Предлагается ввести количество тел для генерации, распределение тел в пространстве (минимальная координата и максимальная по значению координата для x, y и z, ввод через пробел), распределение скоростей тел (минимальное и максимальное значение вектора скорости для x, y и z, ввод через пробел), распределение масс тел (минимальное и максимальное значение массы, ввод через пробел), распределение радиусов тел (минимальное и максимальное значение радиуса, ввод через пробел).
Далее программа выведет число сгенерированных тел (тело может не сгенерироваться, если координаты вместе с радиусом создают коллизию для другого тела) и предлагает выбрать имя файла, который будет сохранен в каталог.
___
### 3 пункт:
Предлагается ввести T, и вероятность разделения тела, после чего выводится время расчета и количество тел для много- и одно-поточной реализации.
___
### 4 пункт:
Сравнение полученных результатов. Если между данными найдется несоответствие для введенного EPSILON, программа выведет номер тела, их данные и разницу между ними. Далее вводя названия для выходных файлов, сохраняются данные расчета для много- и одно-поточной реализации.
___
## Примечание:
> Координаты выражены в _м_, скорости в _м/с_, масса в _кг_, радиусы в _м_. Гравитационная постоянная составляет _6,7*10-3_ для более высокой вероятности соударения тел и демонстрации работы.
