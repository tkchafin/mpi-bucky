#ifndef RAND_H_
#define RAND_H_


class Rand {
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2000, 2003  The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* A version of Marsaglia-MultiCarry */

 public:

  Rand(unsigned int seed1=1234, unsigned int seed2=5678) : I1(seed1), I2(seed2) {}

  void setSeed(unsigned int i1=1234, unsigned int i2=5678) { I1 = i1; I2 = i2; }

  void getSeed(unsigned int& i1, unsigned int& i2) { i1 = I1; i2 = I2; }

  double runif() {
    // Returns a pseudo-random number between 0 and 1.

    I1= 36969*(I1 & 0177777) + (I1>>16);
    I2= 18000*(I2 & 0177777) + (I2>>16);
    return ((I1 << 16)^(I2 & 0177777)) * 2.328306437080797e-10; /* in [0,1) */
  }
 private:
  unsigned int I1, I2;
};

#endif /* RAND_H_ */
