#ifndef CCURVE_H
#define CCURVE_H

#include "../math/vector2d.h"
#include <list>
#include "geometry.h"
#include "cgeom.h"
#include "cubecurve.h"


/**
@brief Квадратная кривая Безье
@author Феклушкин Денис
@date 14 Nov 2008

Квадратная кривая Безье является отрезком обыкновенной квадратной параболы (наклонной).
Задаётся квадратная кривая Безье двумя параметрами: направлением и контрольной точкой.
Внутри класса парабола задана параметрически: с изменением t изменяются x и y.
Кривой Безье при этом принадлежит отрезок этой параболы t [0,1]

Ненужные подробности:
  Любая точка Pt на кривой Безье второго порядка вычисляется по формуле (1):
  Pt = S*(1-t)^2 + 2*C*(1-t)*t + E*t^2
  где:
  t (time) — time-итератор точки;
  S (start) — начальная опорная (узловая) точка (t=0) (anchor point);
  С (control) — управляющая точка (direction point);
  E (end) — конечная опорная (узловая) точка (t=1) (anchor point).

*/

class CCurve : public CGeom
{
public:
    Vector2D direction; ///< Направление кривой (direction vector)
    Vector2D curvature; ///< Кривизна кривой (control vector)

    CCurve();
    CCurve( const CCurve& src );
    CCurve( const Vector2D& direction );
    CCurve( const Vector2D& direction, const Vector2D& curvature );

    ~CCurve();

    /// Подсчёт площади фигуры, ограниченной данной кривой
    Scalar CalculateArea() const;

	/// Подсчитать емкость момента инерции ( чтобы получить момент нужно домножить на плотность )
	Scalar CalcMomentOfInertiaCapacity( Vector2D* axis = 0 ) const;

    /// Аппорксимирует кривую в набор точек
    VectorsList* ConvertToVectors( unsigned num_steps ) const;

    /// Устанавливает цвет отображения данной кривой
    inline void setColor (const Color color) { CCurve::color = color; }

    /// Найти точку, принадлежащую кривой, при заданном t
    Vector2D GetParabolaPoint( Scalar t ) const;

    /// Найти точку касания, соответствующую заданной касательной к параболе
    Scalar GetPointOfOsculation( Vector2D& tangent ) const;

    /// Найти касательную к параболе в заданной точке t
    Vector2D GetTangent( Scalar t ) const;

    /// Найти центр масс фигуры, образованной кривой
    Vector2D CalculateCenterMass() const;

    /// Найти ось параболы
    Vector2D GetParabolaAxis() const;

    /// Найти "габаритный ящик" кривой
    void GetDimensionsBox(Vector2D& v1, Vector2D& v2) const;

    /// Вернуть контрольную точку
    inline Vector2D GetCurvature() const { return curvature; }

    /// Вернуть конец кривой (направление)
    inline Vector2D GetDirection() const { return direction; }

    /// Попадает ли точка в кривую?
    bool IsPointInsideCurve( Vector2D& point ) const;

    /// Найти длину кривой
    Scalar GetCurveLength() const;

    /// Деление кривой на две в заданной точке
    void Split( Scalar t, CCurve& cur1, CCurve& cur2 ) const;

    /// Вогнутая ли кривая?
    bool IsConcave() const;

    /// Поиск пересечения двух парабол
    int Intersection( const CCurve& remote_curve, const Vector2D& remote_coords, Scalar results[4] ) const;

    /// Поворот параболы
    void Rotate( Scalar angle );

protected:
    void Init(); ///< Для внутреннего использования
};

#endif
