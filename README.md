# Атмосфера Земли верхняя ГОСТ Р 25645.166-2004

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![en](https://img.shields.io/badge/lang-en-green.svg)](README.EN.md) [![ru](https://img.shields.io/badge/lang-ru-green.svg)](README.md)

[English version](README.EN.md)

![В разработке](under_construction.png)

*(В доработке, ещё не протестировано!)* Функция для расчёта плотности верхней атмосферы Земли по модели ГОСТ Р 25645.166-2004 для C++, Python и Matlab. Модель ГОСТ Р 25645.166-2004 - российский аналог моделей верхней атмосферы NRLMSISE-00, Яккиа-Боумана JB-2008 и других.

## Use in C++ projects

1) Поместите `atmosGOST_R_25645_166_2004.h` и `atmosGOST_R_25645_166_2004.cpp` в каталог вашего проекта.
2) `#include "atmosGOST_R_25645_166_2004.h"`
3) Плотность рассчитывается вызовом функции `atmosGOST_R_25645_166_2004()` с соответствующими параметрами (см. в комментарии к параметрам функции либо в сам ГОСТ Р 25645.166-2004).

## TODO list

- [ ] Обёртка для Python and Matlab.
- [ ] Оптимизация.
- [ ] Документация.
- [ ] Многоязычная версия репозитория.
- [ ] Адаптация под CMake.
- [ ] Модульное тестирование.
