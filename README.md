# Smirnov-Grubbs Outlier Rejection Technique (Monte Carlo Simulation in MATLAB)

This repository contains a parallel computing MATLAB script designed for large-scale Monte Carlo simulation (10^10 iterations) to calculate the refined and extended critical values for the two-sided Smirnov–Grubbs outlier rejection test (n < 10,000,000).

The repository serves as a supplementary material for the scientific paper:
> **Salmin I.A.** *Improving the Statistical Outlier Rejection Technique for Data Analysis.* (2026).

---

## 🇷🇺 Описание проекта (Russian Description)

В данном репозитории размещен программный скрипт на языке MATLAB, реализующий алгоритм параллельных вычислений методом Монте-Карло (10^10 симуляций). Скрипт предназначен для вычисления уточненных и расширенных таблиц критических значений двустороннего критерия Смирнова–Граббса при объемах выборок n до 10 000 000.

Репозиторий является дополнительным материалом к статье:
> **Салмин И.А.** *Совершенствование методики статистической отбраковки аномальных результатов испытаний.* (2026).

---

## 🛠️ Requirements / Требования

* MATLAB R2021a or newer (R2023–R2026 recommended)
* Parallel Computing Toolbox (for `parfor` multi-core execution)
* Minimum 16 GB RAM (128 GB recommended for large n and k)

---

## 🚀 How to Run / Инструкция по запуску

1. Clone or download this repository.
2. Open `Grubbs_Simulation.m` in MATLAB.
3. Configure simulation parameters if needed:
   * `alphaVal`: Significance level (default is `0.05 * 0.5` for two-sided 5% test).
   * `sampleSizes`: Target sample sizes array (e.g., `3:150`).
   * `numSubSimulations` & `numIteration`: Controls the total Monte Carlo iterations count.
4. Run the script. The results will be displayed in the Command Window and exported automatically to `Grubbs_Critical_Values.xlsx`.

---

## 📑 Citation / Цитирование

If you find this code useful for your research or engineering tasks, please cite the original paper:

```bibtex
@article{salmin2026grubbs,
  author    = {Salmin, I. A.},
  title     = {Improving the Statistical Outlier Rejection Technique for Data Analysis},
  journal   = {Soil Mechanics and Foundation Engineering},
  year      = {2026},
  number    = {12},
  language  = {russian}
}
```

---

## 📄 License / Лицензия

This project is licensed under the **MIT License** - see the LICENSE file for details.
