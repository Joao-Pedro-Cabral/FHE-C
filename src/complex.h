
#ifndef __COMPLEX_H__
#define __COMPLEX_H__

// Real part multiplication
static inline double rel_mul(double a_rel, double a_img, double b_rel,
                             double b_img) {
  return a_rel * b_rel - a_img * b_img;
}

// Imaginary part multiplication
static inline double img_mul(double a_rel, double a_img, double b_rel,
                             double b_img) {
  return a_rel * b_img + a_img * b_rel;
}

#endif  // __COMPLEX_H__
