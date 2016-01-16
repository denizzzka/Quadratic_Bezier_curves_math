#include "ccurve.h"


CCurve::CCurve()
{
    Init();
}


CCurve::CCurve( const CCurve& src )
{
    direction = src.direction;
    curvature = src.curvature;
    Init();
}


CCurve::CCurve( const Vector2D& direction )
{
    CCurve::direction = direction;
    Init();
}


CCurve::CCurve( const Vector2D& direction, const Vector2D& curvature )
{
    CCurve::direction = direction;
    CCurve::curvature = curvature;
    Init();
}


CCurve::~CCurve()
{
}


/*!
    \fn CCurve::Init()
 */
void CCurve::Init()
{
    DisplayDots = true;
    color = ColorCurve;
}


/**
Подсчёт площади фигуры, ограниченной данной кривой и прямой, соединяющей начало и конец кривой.
Площадь кривой положительна если контрольный вектор лежит справа по направлению кривой.
(Аналогично вычислению ориентированной площади многоугольника.)
*/
Scalar CCurve::CalculateArea() const
{
    /*
    Площадь сегмента параболы равна 2/3 от площади треугольника.
    Площадь треугольника равна половине векторного произведения
    двух сторон.
    */
    return (2.0/3.0) * Vector2D::ProductVector(curvature, direction) * 0.5;
}


void CCurve::display() const
{

    // получим "векторизацию" кривой
    unsigned num_steps = static_cast<unsigned>( (direction.length() + curvature.length()) * 3 ) + 2; // FIXME округление может сделать?
    VectorsList* line = ConvertToVectors( num_steps );

    VectorsList::iterator i = line->begin();
    VectorsList::iterator end2 = line->end();

    glColor(color);
    glLineWidth( 2.0f );
    glBegin( GL_LINE_STRIP ); // ломанная линия

    glVertex2f( 0.0, 0.0 ); // первая линия

    while ( i != end2 ) {
      glVertex2f( i->val.x, i->val.y );
      i++;
    }

    glEnd();
    glLineWidth ( 1.0 );

    // проставим ограничительные точки в начале и конце ломанной
    if( DisplayDots ) {

        glPointSize( 4.0 );
        glColor(ColorCurveDots);

        glBegin(GL_POINTS);
          glVertex2f( 0.0, 0.0 );
          glVertex2f( direction.val.x, direction.val.y );
        glEnd();

        glPointSize( 1.0 );

    }

    delete line;

}


/**
   @brief Поиск векторов, приблизительно составляющих кривую.

   Результирующий список может быть меньше количества заказанных шагов!

   @param num_steps максимальное количество шагов минус 2 (начало и конец)
   @return список точек кривой, включая начало и конец
 */
VectorsList* CCurve::ConvertToVectors( unsigned num_steps ) const
{

  // имеет смысл вызывать, чтобы поделить кривую хотя бы на 2 вектора
  assert( num_steps >= 2 );

  // создадим лист
  VectorsList* result_list = new VectorsList;

  result_list->push_back( Vector2D() ); // первый сегмент, точно смотрящий в начало

  // проверим, вдруг наша кривая совсем даже не кривая - сэкономим время
  if ( !curvature.IsNull() ) {

    Scalar t_step = 1.0 / num_steps;
    Scalar t = t_step; // первая точка на кривой, отстоит на один шаг от начала

    while( t < 1.0 ) {

        Vector2D new_point = GetParabolaPoint( t );
        result_list->push_back( new_point );

        t += t_step;
    }

  }

  result_list->push_back( direction ); // последний сегмент, точно смотрящий в конец кривой
  return result_list;

}


/**
 * @param t Коэффициент, "проезжающий" направление от начала до конца на участке [0...1]
 * @return Точка на кривой
 */
Vector2D CCurve::GetParabolaPoint( Scalar t ) const
{
    Vector2D res;

    // формула точки на квадратной кривой Безье: 2 * Control * (1-t)*t + End * t^2
    // она же - параметрически заданная наклонная парабола
    res.val.x = 2.0 * curvature.val.x * (1.0 - t) * t + direction.val.x * pow(t, 2.0);
    res.val.y = 2.0 * curvature.val.y * (1.0 - t) * t + direction.val.y * pow(t, 2.0);

    return res;
}


/*!
    \fn CCurve::GetParabolaExtremum()
 */
Vector2D CCurve::GetParabolaExtremum() const
{
    return GetParabolaPoint( 0.5 ); //TODO: неправильная точка, исправить
}


/**
@param tangent Касательная к параболе
@return Точка t, соответствующая касательной
*/
Scalar CCurve::GetPointOfOsculation( Vector2D& tangent ) const
{
  // при нулевой длине кривой и нулевой кривизне точка касания равна нулю
  if( direction.IsNull() && curvature.IsNull() ) return 0.0;

  Scalar k;
  if( tangent.val.x ) k = tangent.val.y / tangent.val.x; else k = 0.0; // x==0, поэтому k=0

  // k = y'(t) / x'(t), отсюда t:
  Scalar t = (2.0 * curvature.val.y - 2.0 * curvature.val.x * k) /
              (4.0 * curvature.val.y - 2.0 * direction.val.y - 4.0 * k * curvature.val.x +
               2.0 * direction.val.x * k );
  return t;
}


Vector2D CCurve::GetTangent( Scalar t ) const
{
  Vector2D res;

  // если кривая совсем не кривая, то, избегая возвращения ноля при t=0, вернём просто
  // параллельный вектор
  if( curvature.IsNull() ) {
      res = direction;
  } else {
      // это просто производная формул кривой Безье из CCurve::GetParabolaPoint()
      res.val.x = 2.0 * curvature.val.x - 4.0 * curvature.val.x * t + 2.0 * direction.val.x * t;
      res.val.y = 2.0 * curvature.val.y - 4.0 * curvature.val.y * t + 2.0 * direction.val.y * t;
  }

  assert( !res.IsNull() );
  return res;
}


Vector2D CCurve::CalculateCenterMass() const
{
    // середина направления
    Vector2D dir_half = 0.5 * direction;

    // TODO считаем на глазок, без формул, похоже на правду но не доказано :(
    // вернём треть от высоты параболы
    return dir_half + GetParabolaAxis() * (1.0/3.0);
}


/*!
    \fn CCurve::GetParabolaAxis()
 */
Vector2D CCurve::GetParabolaAxis() const
{
    return GetParabolaExtremum() - direction * 0.5;
}


/*!
    \fn CCurve::GetDimensionsBox(Vector2D& v1, Vector2D& v2)
 */
void CCurve::GetDimensionsBox(Vector2D& v1, Vector2D& v2) const
{
  // для начала закинем начало
  Vector2D start;
  Vector2D::MergeVector ( start, &v1, &v2 );

  // конец
  Vector2D::MergeVector ( direction, &v1, &v2 );

  // половину кривизны
  Vector2D c( curvature * 0.5 );
  Vector2D::MergeVector ( c, &v1, &v2 );

  // и половину другой стороны кривизны
  Vector2D b( (curvature - direction) * 0.5 + direction );
  Vector2D::MergeVector ( b, &v1, &v2 );
}


bool CCurve::IsPointInsideCurve( Vector2D& point ) const
{
    // направление на точку из вершины curvature
    Vector2D to_point = point - curvature;

    // точка пересечения с direction;
    Vector2D cross = direction.CrossingPoint( &to_point, &curvature );

    // найдём соответствующую этой точке t
    Scalar t = cross / direction;

    // если t вышла за пределы то точка уже точно не попадает
    if( t < 0 || t > 1 ) return false;

    // найдём соответствующую t точку на кривой относительно cross
    Vector2D etalon = GetParabolaPoint( t ) - cross;

    // найдём вектор от cross к изучаемой точке
    Vector2D dir2point = point - cross;

    // если точка лежит на direction то она принадлежит площади кривой
    if( dir2point.IsNull() ) return true;

    // проверим, в одну ли сторону смотрят векторы
    if( !dir2point.SectorCompare( &etalon ) ) return false;

    // проверяем длины векторов, делаем выводы :)
    if( dir2point.length() > etalon.length() ) return false; else return true;
}


/**
 @return длина кривой
 */
Scalar CCurve::GetCurveLength() const
{
    Vector2D control2end = direction - curvature;
    Scalar vec_len = direction.length() + curvature.length() + control2end.length();
    return 0.5 * vec_len;
}



/**
Оригинальная кривая остаётся нетронутой.

@param t Точка, в которой будет разделена кривая
@param cur1 Сюда будет помещена первая по направлению полученная кривая
@param cur2 Сюда будет помещена вторая по направлению полученная кривая
*/
void CCurve::Split( Scalar t, CCurve& cur1, CCurve& cur2 ) const
{
    assert( t >= 0 );
    assert( t <= 1 );

    // координаты точки разделения
    cur1.direction = GetParabolaPoint( t );

    // касательная в точке разделения
    Vector2D tangent = GetTangent( t );

    // контролька первой кривой
    cur1.curvature = curvature.CrossingPoint( &tangent, &cur1.direction );

    // направление второй кривой
    cur2.direction = direction - cur1.direction;

    // недостающая сторона треугольника
    Vector2D c2 = direction - curvature;

    // вторая контролька
    cur2.curvature = tangent.CrossingPoint( &c2, &cur2.direction );
}


/**
 @brief Проверка кривой на вогнутость.

 @return если крива я вогнутая - true, прямая или выпуклая - false
*/
bool CCurve::IsConcave() const
{
    return direction.Classify( curvature ) == Vector2D::LEFT;
}

/**
 @param axis Необязательный параметр - ось вращения
 @author Александр Оплеснин <gasparfx@gmail.com>
*/
Scalar CCurve::CalcMomentOfInertiaCapacity( Vector2D* axis ) const
{
	// Если объект слишком мал - не считаем
	Scalar area = fabs( CalculateArea() );
	if ( area < 0.0000001 ){
		return 0;
	}
	// Сюда будем результат запоминать
	Scalar momInerCap = 0;

	// Ось вращения
	Vector2D axisOfRotation;
	// Центр масс фигуры
	Vector2D centrMass = CCurve::CalculateCenterMass();

	// Если ось не была задана - берем центр масс
	if ( !axis ){
		axisOfRotation = centrMass;
	} else {
		axisOfRotation = *axis;
	}
        Scalar dirAngle = atan2(direction.val.y, direction.val.x); ///< Угол между direction и Ox(против часовой)
        Scalar curAngle = atan2(curvature.val.y, curvature.val.x); ///< Угол между curvature и Ox(против часовой)
        Scalar newCurAngle = curAngle - dirAngle;
        Vector2D d(direction.length(), 0);
        Vector2D c(cos( newCurAngle ) * curvature.length(), sin( newCurAngle ) * curvature.length() );

        Scalar Cx = c.val.x;
        Scalar Cy = c.val.y;
        Scalar Dx = d.val.x;

		// Емкость момента инерции относительно начала
		momInerCap = Cy*Dx*(4*Cx*Cx+10*Dx*Cx+4*Cy*Cy+15*Dx*Dx)/210.0; /// TODO Потестить формулу

		// Емкость момента инерции относительно центра масс
		momInerCap = fabs( momInerCap ) - ( area * centrMass.length() * centrMass.length() );


		momInerCap = momInerCap + ( area * ( centrMass - axisOfRotation ).length()
										 * ( centrMass - axisOfRotation ).length() );

		return momInerCap;
}



/**
 @brief Поиск пересечений двух парабол

 @param remote_curve - парабола, с которой проверяется пересечение
 @param remote_coords - координаты второй параболы
 @param results - массив точек пересечения, заданных t
         (t идут в сторону уменьшения)
 @return количество пересечений (0, 1, 2 или 4)
         или -1 - кривые совпадают
*/
int CCurve::Intersection( const CCurve& remote_curve, const Vector2D& remote_coords, Scalar results[4] ) const
{
    const CCurve &c2 = remote_curve;

    const Vector2D S = -1 * remote_coords;

    Vector2D D1 = direction - 2.0 * curvature;
    Vector2D D2 = c2.direction - 2.0 * c2.curvature;

    Scalar Qc1 = 2.0 * curvature.val.x * D2.val.y - 2.0 * curvature.val.y * D2.val.x;
    Scalar Qc2 = 2.0 * c2.curvature.val.x * D2.val.y - 2.0 * c2.curvature.val.y * D2.val.x;

    Scalar Qd = D1.val.x * D2.val.y - D1.val.y * D2.val.x;
    Scalar Qs = S.val.x * D2.val.y - S.val.y * D2.val.x;

    Scalar Rc = Qc1 / Qc2;
    Scalar Rd = Qd / Qc2;
    Scalar Rs = Qs / Qc2;

    // Коэффициенты для вызова функции решения уравнения 4 степени поместим в массив:
    Scalar abcde[5] = {
        Rd * Rd, //a
        2 * Rd * Rc, //b
        (Rc*Rc + 2.0*Rd*Rs - 2.0*(c2.curvature.val.x*D1.val.y - c2.curvature.val.y * D1.val.x)/Qc2), //c
        (2.0*Rs*Rc + 4.0*(curvature.val.x * c2.curvature.val.y - c2.curvature.val.x * curvature.val.y)/Qc2), //d
        (Rs*Rs + 2.0*((S.val.x*c2.curvature.val.y - S.val.y*c2.curvature.val.x)/Qc2)) //e
    };

    int res = fast_math::QuarticRR( abcde, results );

    // определим количество полученных точек
    switch( res ){
        case 6: // 4 реальных корня
            return 4;
            break;
        case 2: // 2 реальных корня
            return 2;
            break;
        case -1: // 1 корень
            return 1;
            break;
        case -2: // 0 корней
            return 0;
            break;
        case -3: // бесконечное число корней
            return -1;
        default:
            PR( res );
            PR( remote_coords );
            PR( direction );
            PR( curvature);
            PR( c2.direction );
            PR( c2.curvature );
            assert(false);
            break;
    }

    return -8; // никогда  используется
}


/**
 @param angle - поворота параболы в радианах
        (положительное значение соответствует повороту по часовой стрелке)
*/
void CCurve::Rotate( Scalar angle )
{
    direction.Rotate( angle );
    curvature.Rotate( angle );
}
