#ifndef TEST_BOUNDARIES
#define TEST_BOUNDARIES
//#############################################################################//
// Parametric Boundaries
//#############################################################################//
// Parametric function for boundary part 1

/// \short Specialisation of CurvilineGeomObject for half a circle.
class CurvilineQuadratic : public CurvilineGeomObject
{
public:
  /// Constructor
  CurvilineQuadratic() {};

  /// Destructor
  ~CurvilineQuadratic() {};

 /// Broken copy constructor
 CurvilineQuadratic(const CurvilineQuadratic& dummy) 
  { 
   BrokenCopy::broken_copy("CurvilineQuadratic");
  } 
 
 /// Broken assignment operator
 void operator=(const CurvilineQuadratic&) 
  {
   BrokenCopy::broken_assign("CurvilineQuadratic");
  }
 /// Parametric function
 void position(const Vector<double>& s, Vector<double>& x) const
 { x[0] =-s[0];  x[1] =-s[0]*s[0]+0.75;};
 
 /// Derivative of parametric function
 void dposition(const Vector<double>& s, Vector<double>& dx) const
 { dx[0] =-1;  dx[1] =-2*s[0];};
 
 /// Derivative of parametric function
 void d2position(const Vector<double>& s, Vector<double>& dx) const 
 { dx[0] =0;  dx[1] =-2;};
};


//#############################################################################//
// Parametric Boundaries
//#############################################################################//
// Parametric function for boundary part 1

/// \short Specialisation of CurvilineGeomObject for half a circle.
class CurvilineCubic : public CurvilineGeomObject
{
public:
  /// Constructor
  CurvilineCubic() {};

  /// Destructor
  ~CurvilineCubic() {};

 /// Broken copy constructor
 CurvilineCubic(const CurvilineCubic& dummy) 
  { 
   BrokenCopy::broken_copy("CurvilineCubic");
  } 
 
 /// Broken assignment operator
 void operator=(const CurvilineCubic&) 
  {
   BrokenCopy::broken_assign("CurvilineCubic");
  }
 /// Parametric function
 void position(const Vector<double>& s, Vector<double>& x) const
 { x[0] = s[0];  x[1] =s[0]*(s[0]*s[0]-1/4.)+sqrt(3)/2.;;};
 
 /// Derivative of parametric function
 void dposition(const Vector<double>& s, Vector<double>& dx) const
 { dx[0] = 1;  dx[1] =(3.0*s[0]*s[0]-1/4.);};
 
 /// Derivative of parametric function
 void d2position(const Vector<double>& s, Vector<double>& dx) const 
 { dx[0] =0;  dx[1] =(6.0*s[0]);};
};

/// \short Specialisation of CurvilineGeomObject for half a circle.
class CurvilineCircleRight : public CurvilineGeomObject
{
public:
 /// \short Constructor: Pass dimension of geometric object (# of Eulerian
 /// coords = # of Lagrangian coords; no time history available/needed)
 CurvilineCircleRight() : CurvilineGeomObject()
  { }

 /// \short Constructor: pass # of Eulerian and Lagrangian coordinates
 /// and pointer to time-stepper which is used to handle the
 /// position at previous timesteps and allows the evaluation
 /// of veloc/acceleration etc. in cases where the GeomData
 /// varies with time.
 CurvilineCircleRight(TimeStepper* time_stepper_pt) : CurvilineGeomObject(time_stepper_pt)
   { }

 /// Broken copy constructor
 CurvilineCircleRight(const CurvilineCircleRight& dummy) 
  { 
   BrokenCopy::broken_copy("CurvilineCircleRight");
  } 
 
 /// Broken assignment operator
 void operator=(const CurvilineCircleRight&) 
  {
   BrokenCopy::broken_assign("CurvilineCircleRight");
  }

 /// (Empty) destructor
 virtual ~CurvilineCircleRight(){}

 /// \short Position Vector w.r.t. to zeta: 
 virtual void position(const Vector<double>& zeta, 
                        Vector<double> &r) const
  { r[0] = std::cos(zeta[0]);  r[1] = std::sin(zeta[0]);}

 /// \short Derivative of position Vector w.r.t. to zeta: 
 virtual void dposition(const Vector<double>& zeta, 
                        Vector<double> &drdzeta) const
  { drdzeta[0] =-std::sin(zeta[0]);  drdzeta[1] = std::cos(zeta[0]);}


 /// \short 2nd derivative of position Vector w.r.t. to coordinates: 
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ = 
 /// ddrdzeta(alpha,beta,i). 
 /// Evaluated at current time.
 virtual void d2position(const Vector<double>& zeta, 
                         Vector<double> &drdzeta) const
  { drdzeta[0] =-std::cos(zeta[0]);  drdzeta[1] =-std::sin(zeta[0]);}

 /// Get s from x for part 0 of the boundary (inverse mapping - for convenience)
 double get_zeta(const Vector<double>& x)
 {
 // The arc length (parametric parameter) for the upper semi circular arc
  return atan2(x[0],-x[1]);
 }
  

};

/// \short Specialisation of CurvilineGeomObject for half a circle.
class CurvilineStraightLine : public CurvilineGeomObject
{
public:
 /// \short Constructor: Pass dimension of geometric object (# of Eulerian
 /// coords = # of Lagrangian coords; no time history available/needed)
 CurvilineStraightLine() : CurvilineGeomObject(), Mx(0), Cx(0), My(0), Cy(0)
  {}

 /// \short Constructor: Pass dimension of geometric object (# of Eulerian
 /// coords = # of Lagrangian coords; no time history available/needed)
 CurvilineStraightLine(const Vector<double>& xcoeffs, const Vector<double>& ycoeffs) 
   : CurvilineGeomObject(), Mx(xcoeffs[1]), Cx(xcoeffs[0]), My(ycoeffs[1]), Cy(ycoeffs[0])
  {}

 /// \short Constructor: pass # of Eulerian and Lagrangian coordinates
 /// and pointer to time-stepper which is used to handle the
 /// position at previous timesteps and allows the evaluation
 /// of veloc/acceleration etc. in cases where the GeomData
 /// varies with time.
 CurvilineStraightLine(TimeStepper* time_stepper_pt) : CurvilineGeomObject(time_stepper_pt)
   { }

 /// Broken copy constructor
 CurvilineStraightLine(const CurvilineStraightLine& dummy) 
  { 
   BrokenCopy::broken_copy("CurvilineStraightLine");
  } 
 
 /// Broken assignment operator
 void operator=(const CurvilineStraightLine&) 
  {
   BrokenCopy::broken_assign("CurvilineStraightLine");
  }

 /// (Empty) destructor
 virtual ~CurvilineStraightLine(){}

 /// \short Position Vector w.r.t. to zeta: 
 virtual void position(const Vector<double>& zeta, 
                        Vector<double> &r) const
  { r[0] = Mx * zeta[0] + Cx;  r[1] = My * zeta[0] + Cy;}

 /// \short Derivative of position Vector w.r.t. to zeta: 
 virtual void dposition(const Vector<double>& zeta, 
                        Vector<double> &drdzeta) const
  { drdzeta[0] = Mx;  drdzeta[1] = My;}

 /// \short 2nd derivative of position Vector w.r.t. to coordinates: 
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ = 
 /// ddrdzeta(alpha,beta,i). 
 /// Evaluated at current time.
 virtual void d2position(const Vector<double>& zeta, 
                         Vector<double> &drdzeta) const
  { drdzeta[0] = 0;  drdzeta[1] =0;}

 /// Get s from x for part 0 of the boundary (inverse mapping - for convenience)
 double get_zeta(const Vector<double>& x)
 {
 // The arc length (parametric parameter) for the upper semi circular arc
  return (x[0]-Cx) / Mx;
 }
 
private:
 /// The Gradient of the x component wrt to parametric coordinate
 double Mx; 
 /// The offset of the x component
 double Cx; 
 /// The Gradient of the y component wrt to parametric coordinate
 double My; 
 /// The offset of the x component
 double Cy; 
};
#endif
