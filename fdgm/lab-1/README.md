# Лабораторная работа №1 — метод конечных разностей для ОДУ второго порядка

Краевая задача третьего рода

    -u'' + p(x)u' + q(x)u = f(x),   x ∈ [a,b]
    -α_a u'(a) + β_a u(a) = γ_a
     α_b u'(b) + β_b u(b) = γ_b

Данные: `[a,b] = [0,2]`, `p = 2 + cos πx`, `q = 1 + x²`,
`α_a = β_a = α_b = 1`, `β_b = 2`, точное решение `u* = x sin 2πx + e^{-x}`.

## Файлы

    code/problem.py    постановка: p, q, u*, f*, γ*
    code/schemes.py    сборка СЛАУ (схемы O(h) и O(h²)) и прогонка
    code/run.py        расчеты, таблицы, графики
    report/report.tex  отчет (pdfLaTeX)
    report/figures/    графики

## Запуск

    python code/run.py

Печатает таблицы в консоль и пересобирает `report/figures/*.pdf`.
Нужны `numpy` и `matplotlib`.

Отчет компилируется в Overleaf: загрузить папку `report/` целиком.
