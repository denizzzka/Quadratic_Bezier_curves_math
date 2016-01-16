(RU) Квадратная кривая Безье (кривая Безье второго порядка)
------

Квадратная кривая Безье является отрезком обыкновенной (наклонной) квадратной параболы.
Задаётся квадратная кривая Безье двумя параметрами: направлением и контрольной точкой.
Внутри класса парабола задана параметрически: с изменением t изменяются x и y.
Кривой Безье при этом принадлежит отрезок этой параболы t [0,1]

<pre>
Любая точка Pt на кривой Безье второго порядка вычисляется по формуле:
Pt = S*(1-t)^2 + 2*C*(1-t)*t + E*t^2
где:
t (time) — time-итератор точки;
S (start) — начальная опорная (узловая) точка (t=0) (anchor point);
С (control) — управляющая точка (direction point);
E (end) — конечная опорная (узловая) точка (t=1) (anchor point).
</pre>
Этот код может делать:

* Подсчёт площади фигуры, ограниченной данной кривой
* Подсчитать емкость момента инерции (чтобы получить момент нужно домножить на плотность)
* Аппорксимирует кривую в набор точек
* Найти точку, принадлежащую кривой, при заданном t
* Найти точку касания, соответствующую заданной вектором касательной к параболе
* Найти касательную к параболе в заданной точке t
* Найти центр масс фигуры, образованной кривой
* Найти ось параболы
* Найти "габаритный ящик" кривой
* Найти длину кривой
* Деление кривой на две кривые в заданной точке
* Поиск пересечения двух парабол


(EN) Quadratic Bezier curve
====

This code can do:

* Calculation area of ​​the figure bounded by a given curve
* Calculate the capacity of the moment of inertia (to get moment of inertia you need to multiply it by density)
* Curve approximation into a set of points
* Find a point belonging to a curve at a given t
* Find a point of contact corresponding to a given tangent vector to the parabola
* Find the tangent to the parabola at the given point t
* Find the center of mass of the figure formed by the curve
* Find the axis of the parabola
* Find "bounding box"
* Find the length of the curve
* The division of the curve at the predetermined point into two curves
* Search the intersection of two parabolas
