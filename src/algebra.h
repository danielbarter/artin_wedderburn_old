typedef struct algebra {

  // we use the notation e_1, e_2,..., e_dimension for the basis vectors
  int dimension;


  // unit of algebra. It satisfies x unit = unit x = x for all x
  // the unit has length dimension
  double *unit;


  // multiplication tensor of the algebra. This is a 3-tensor.
  // multiplication tensor has length dimension along each axis
  // multiplication[i][j][k] is the coefficient of e_k in e_i e_j
  // to minimize cache misses, we don't store arrays of intermediate pointers
  // instead we compute offsets on the fly
  // this means you can't use the syntax multiplication[i][j][k]
  // instead, use the method get_multiplication_3

  // the multiplication values are stored in row major order
  // as you step through adjacent memory locations
  // the first index varies most slowly
  // and the last index varies most quickly
  // this means that left multiplication matrices x -> a x
  // are contiguious blocks in memory

  // if dimension = 100, then multiplcation is 8MB in size
  double *multiplication;


  // true if we have a *-structure
  bool is_star_algebra;

  // star[i][j] is the coefficient of e_j in e_i^*
  // again, to minimize cache misses, we don't store intermediate pointers
  // star values are stored in row major order
  // as you step through adjacent memory locations
  // the first index varies most slowly
  // and the last index varies most quickly
  // if is_star_algebra == false, then star == NULL
  double *star;

} Algebra;

// getters and setters for multiplication[i][j][k]
double get_multiplication_3(Algebra *algebra, int i, int j, int k);
void set_multiplication_3(Algebra *algebra, double v, int i, int j, int k);

// getters and setters for multiplication[i][j]
double *get_multiplication_2(Algebra *algebra, int i, int j);
void set_multiplication_2(Algebra *algebra, double *vp, int i, int j);

// getters and setters for multiplication[i]
// left multiplication matrix x -> a x is returned and set in row major order
double *get_multiplication_1(Algebra *algebra, int i);
void set_multiplication_1(Algebra *algebra, double *vp, int i);
