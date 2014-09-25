/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Christian Schulte <schulte@gecode.org>
 *
 *  Copyright:
 *     Christian Schulte, 2005
 *
 *  Last modified:
 *     $Date: 2013-08-29 16:05:54 +0200 (to, 29 aug 2013) $ by $Author: schulte $
 *     $Revision: 13993 $
 *
 *  This file is part of Gecode, the generic constraint
 *  development environment:
 *     http://www.gecode.org
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include "test/int.hh"

#include <cmath>
#include <algorithm>

#include <gecode/minimodel.hh>

namespace Test { namespace Int {

   /// %Tests for arithmetic constraints
   namespace Arithmetic {

     /**
      * \defgroup TaskTestIntArithmetic Arithmetic constraints
      * \ingroup TaskTestInt
      */
     //@{
     /// %Test for multiplication constraint
     class MultXYZ : public Test {
     public:
       /// Create and register test
       MultXYZ(const std::string& s, const Gecode::IntSet& d,
               Gecode::IntConLevel icl)
         : Test("Arithmetic::Mult::XYZ::"+str(icl)+"::"+s,3,d,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         double d0 = static_cast<double>(x[0]);
         double d1 = static_cast<double>(x[1]);
         double d2 = static_cast<double>(x[2]);
         return d0*d1 == d2;
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::mult(home, x[0], x[1], x[2], icl);
       }
     };

     /// %Test for multiplication constraint with shared variables
     class MultXXY : public Test {
     public:
       /// Create and register test
       MultXXY(const std::string& s, const Gecode::IntSet& d,
               Gecode::IntConLevel icl)
         : Test("Arithmetic::Mult::XXY::"+str(icl)+"::"+s,2,d,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         double d0 = static_cast<double>(x[0]);
         double d1 = static_cast<double>(x[0]);
         double d2 = static_cast<double>(x[1]);
         return d0*d1 == d2;
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::mult(home, x[0], x[0], x[1], icl);
       }
     };

     /// %Test for multiplication constraint with shared variables
     class MultXYX : public Test {
     public:
       /// Create and register test
       MultXYX(const std::string& s, const Gecode::IntSet& d,
               Gecode::IntConLevel icl)
         : Test("Arithmetic::Mult::XYX::"+str(icl)+"::"+s,2,d,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         double d0 = static_cast<double>(x[0]);
         double d1 = static_cast<double>(x[1]);
         double d2 = static_cast<double>(x[0]);
         return d0*d1 == d2;
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::mult(home, x[0], x[1], x[0], icl);
       }
     };

     /// %Test for multiplication constraint with shared variables
     class MultXYY : public Test {
     public:
       /// Create and register test
       MultXYY(const std::string& s, const Gecode::IntSet& d,
               Gecode::IntConLevel icl)
         : Test("Arithmetic::Mult::XYY::"+str(icl)+"::"+s,2,d,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         double d0 = static_cast<double>(x[0]);
         double d1 = static_cast<double>(x[1]);
         double d2 = static_cast<double>(x[1]);
         return d0*d1 == d2;
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::mult(home, x[0], x[1], x[1], icl);
       }
     };

     /// %Test for multiplication constraint with shared variables
     class MultXXX : public Test {
     public:
       /// Create and register test
       MultXXX(const std::string& s, const Gecode::IntSet& d,
               Gecode::IntConLevel icl)
         : Test("Arithmetic::Mult::XXX::"+str(icl)+"::"+s,1,d,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         double d0 = static_cast<double>(x[0]);
         double d1 = static_cast<double>(x[0]);
         double d2 = static_cast<double>(x[0]);
         return d0*d1 == d2;
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::mult(home, x[0], x[0], x[0], icl);
       }
     };

     /// %Test for squaring constraint
     class SqrXY : public Test {
     public:
       /// Create and register test
       SqrXY(const std::string& s, const Gecode::IntSet& d,
             Gecode::IntConLevel icl)
         : Test("Arithmetic::Sqr::XY::"+str(icl)+"::"+s,2,d,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         double d0 = static_cast<double>(x[0]);
         double d1 = static_cast<double>(x[1]);
         return d0*d0 == d1;
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::sqr(home, x[0], x[1], icl);
       }
     };
     
     /// %Test for squaring constraint with shared variables
     class SqrXX : public Test {
     public:
       /// Create and register test
       SqrXX(const std::string& s, const Gecode::IntSet& d,
             Gecode::IntConLevel icl)
         : Test("Arithmetic::Sqr::XX::"+str(icl)+"::"+s,1,d,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         double d0 = static_cast<double>(x[0]);
         return d0*d0 == d0;
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::sqr(home, x[0], x[0], icl);
       }
     };

     /// %Test for square root constraint
     class SqrtXY : public Test {
     public:
       /// Create and register test
       SqrtXY(const std::string& s, const Gecode::IntSet& d,
              Gecode::IntConLevel icl)
         : Test("Arithmetic::Sqrt::XY::"+str(icl)+"::"+s,2,d,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         double d0 = static_cast<double>(x[0]);
         double d1 = static_cast<double>(x[1]);
         return (d0 >= 0) && (d0 >= d1*d1) && (d0 < (d1+1)*(d1+1));
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::sqrt(home, x[0], x[1], icl);
       }
     };

     /// %Test for square root constraint with shared variables
     class SqrtXX : public Test {
     public:
       /// Create and register test
       SqrtXX(const std::string& s, const Gecode::IntSet& d,
              Gecode::IntConLevel icl)
         : Test("Arithmetic::Sqrt::XX::"+str(icl)+"::"+s,1,d,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         double d0 = static_cast<double>(x[0]);
         return (d0 >= 0) && (d0 >= d0*d0) && (d0 < (d0+1)*(d0+1));
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::sqrt(home, x[0], x[0], icl);
       }
     };

     /// %Test for power constraint
     class PowXY : public Test {
     protected:
       /// The exponent
       int n;
     public:
       /// Create and register test
       PowXY(const std::string& s, int n0, const Gecode::IntSet& d,
             Gecode::IntConLevel icl)
         : Test("Arithmetic::Pow::XY::"+str(n0)+"::"+str(icl)+"::"+s,
                2,d,false,icl), n(n0) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         long long int p = 1;
         for (int i=0; i<n; i++) {
           p *= x[0];
           if ((p < Gecode::Int::Limits::min) || 
               (p > Gecode::Int::Limits::max))
             return false;
         }
         return p == x[1];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         using namespace Gecode;
         if (n > 4)
           pow(home, x[0], n, x[1], icl);
         else
           rel(home, expr(home, pow(x[0],n), icl), IRT_EQ, x[1], icl);
       }
     };

     /// %Test for power constraint with shared variables
     class PowXX : public Test {
     protected:
       /// The exponent
       int n;
     public:
       /// Create and register test
       PowXX(const std::string& s, int n0, const Gecode::IntSet& d,
             Gecode::IntConLevel icl)
         : Test("Arithmetic::Pow::XX::"+str(n0)+"::"+str(icl)+"::"+s,
                1,d,false,icl), n(n0) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         long long int p = 1;
         for (int i=0; i<n; i++) {
           p *= x[0];
           if ((p < Gecode::Int::Limits::min) || 
               (p > Gecode::Int::Limits::max))
             return false;
         }
         return p == x[0];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::pow(home, x[0], n, x[0], icl);
       }
     };

     bool powgr(int n, long long int r, int x) {
       assert(r >= 0);
       long long int y = r;
       long long int p = 1;
       do {
         p *= y; n--;
         if (p > x)
           return true;
       } while (n > 0);
       return false;
     }

     int fnroot(int n, int x) {
       if (x < 2)
         return x;
       /*
        * We look for l such that: l^n <= x < (l+1)^n
        */
       long long int l = 1;
       long long int u = x;
       do {
         long long int m = (l + u) >> 1;
         if (powgr(n,m,x)) u=m; else l=m;
       } while (l+1 < u);
       return static_cast<int>(l);
     }

     bool powle(int n, long long int r, int x) {
       assert(r >= 0);
       long long int y = r;
       long long int p = 1;
       do {
         p *= y; n--;
         if (p >= x)
           return false;
       } while (n > 0);
       assert(y < x);
       return true;
     }

     int cnroot(int n, int x) {
       if (x < 2)
         return x;
       /*
        * We look for u such that: (u-1)^n < x <= u^n
        */
       long long int l = 1;
       long long int u = x;
       do {
         long long int m = (l + u) >> 1;
         if (powle(n,m,x)) l=m; else u=m;
       } while (l+1 < u);
       return static_cast<int>(u);
     }

     /// %Test for nroot constraint
     class NrootXY : public Test {
     protected:
       /// The root index
       int n;
       /// Floor 
     public:
       /// Create and register test
       NrootXY(const std::string& s, int n0, const Gecode::IntSet& d,
             Gecode::IntConLevel icl)
         : Test("Arithmetic::Nroot::XY::"+str(n0)+"::"+str(icl)+"::"+s,
                2,d,false,icl), n(n0) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         if (n == 1)
           return x[0] == x[1];
         if ((n % 2 == 0) && ((x[0] < 0) || (x[1] < 0)))
           return false;
         int r = (x[0] < 0) ? -cnroot(n,-x[0]) : fnroot(n,x[0]);
         return r == x[1]; 
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         using namespace Gecode;
         if (n > 4)
           nroot(home, x[0], n, x[1], icl);
         else
           rel(home, expr(home, nroot(x[0],n), icl), IRT_EQ, x[1], icl);
       }
     };

     /// %Test for nroot constraint with shared variables
     class NrootXX : public Test {
     protected:
       /// The root index
       int n;
     public:
       /// Create and register test
       NrootXX(const std::string& s, int n0, const Gecode::IntSet& d,
               Gecode::IntConLevel icl)
         : Test("Arithmetic::Nroot::XX::"+str(n0)+"::"+str(icl)+"::"+s,
                1,d,false,icl), n(n0) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         if (n == 1)
           return true;
         if (n % 2 == 0) {
           return (x[0] >= 0) && (x[0] <= 1);
         } else {
           return (x[0] >= -2) && (x[0] <= 1);
         }
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::nroot(home, x[0], n, x[0], icl);
       }
     };

     /// %Test for division/modulo constraint
     class DivMod : public Test {
     private:
       /// Return the absolute value of \a a
       static int abs(int a) { return a<0 ? -a:a; }
       /// Return the sign of \a a
       static int sgn(int a) { return a<0 ? -1:1; }
     public:
       /// Create and register test
       DivMod(const std::string& s, const Gecode::IntSet& d)
         : Test("Arithmetic::DivMod::"+s,4,d) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return x[0] == x[1]*x[2]+x[3] &&
                abs(x[3]) < abs(x[1]) &&
                (x[3] == 0 || sgn(x[3]) == sgn(x[0]));
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::divmod(home, x[0], x[1], x[2], x[3]);
       }
     };

     /// %Test for division constraint
     class Div : public Test {
     public:
       /// Create and register test
       Div(const std::string& s, const Gecode::IntSet& d)
         : Test("Arithmetic::Div::"+s,3,d) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         if (x[1] == 0)
           return false;
         int divsign = (x[0] / x[1] < 0) ? -1 : 1;
         int divresult =
           divsign *
           static_cast<int>(floor(static_cast<double>(std::abs(x[0]))/
                                  static_cast<double>(std::abs(x[1]))));
         return x[2] == divresult;
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::div(home, x[0], x[1], x[2]);
       }
     };

     /// %Test for modulo constraint
     class Mod : public Test {
     public:
       /// Create and register test
       Mod(const std::string& s, const Gecode::IntSet& d)
         : Test("Arithmetic::Mod::"+s,3,d) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         if (x[1] == 0)
           return false;
         int divsign = (x[0] / x[1] < 0) ? -1 : 1;
         int divresult =
           divsign *
           static_cast<int>(floor(static_cast<double>(std::abs(x[0]))/
                                  static_cast<double>(std::abs(x[1]))));
         return x[0] == x[1]*divresult+x[2];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::mod(home, x[0], x[1], x[2]);
       }
     };

     /// %Test for absolute value constraint
     class AbsXY : public Test {
     public:
       /// Create and register test
       AbsXY(const std::string& s, const Gecode::IntSet& d,
             Gecode::IntConLevel icl)
         : Test("Arithmetic::Abs::XY::"+str(icl)+"::"+s,2,d,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         double d0 = static_cast<double>(x[0]);
         double d1 = static_cast<double>(x[1]);
         return (d0<0 ? -d0 : d0) == d1;
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::abs(home, x[0], x[1], icl);
       }
     };

     /// %Test for absolute value constraint with shared variables
     class AbsXX : public Test {
     public:
       /// Create and register test
       AbsXX(const std::string& s, const Gecode::IntSet& d,
             Gecode::IntConLevel icl)
         : Test("Arithmetic::Abs::XX::"+str(icl)+"::"+s,1,d,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         double d0 = static_cast<double>(x[0]);
         double d1 = static_cast<double>(x[0]);
         return (d0<0 ? -d0 : d0) == d1;
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::abs(home, x[0], x[0], icl);
       }
     };

     /// %Test for binary minimum constraint
     class MinXYZ : public Test {
     public:
       /// Create and register test
       MinXYZ(const std::string& s, const Gecode::IntSet& d,
              Gecode::IntConLevel icl)
         : Test("Arithmetic::Min::Bin::XYZ::"+str(icl)+"::"+s,3,d,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::min(x[0],x[1]) == x[2];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::min(home, x[0], x[1], x[2], icl);
       }
     };

     /// %Test for binary minimum constraint with shared variables
     class MinXXY : public Test {
     public:
       /// Create and register test
       MinXXY(const std::string& s, const Gecode::IntSet& d,
              Gecode::IntConLevel icl)
         : Test("Arithmetic::Min::Bin::XYX::"+str(icl)+"::"+s,2,d) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::min(x[0],x[0]) == x[1];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::min(home, x[0], x[0], x[1], icl);
       }
     };

     /// %Test for binary minimum constraint with shared variables
     class MinXYX : public Test {
     public:
       /// Create and register test
       MinXYX(const std::string& s, const Gecode::IntSet& d,
              Gecode::IntConLevel icl)
         : Test("Arithmetic::Min::Bin::XYX::"+str(icl)+"::"+s,2,d) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::min(x[0],x[1]) == x[0];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::min(home, x[0], x[1], x[0], icl);
       }
     };

     /// %Test for binary minimum constraint with shared variables
     class MinXYY : public Test {
     public:
       /// Create and register test
       MinXYY(const std::string& s, const Gecode::IntSet& d,
              Gecode::IntConLevel icl)
         : Test("Arithmetic::Min::Bin::XYY::"+str(icl)+"::"+s,2,d) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::min(x[0],x[1]) == x[1];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::min(home, x[0], x[1], x[1], icl);
       }
     };

     /// %Test for binary minimum constraint with shared variables
     class MinXXX : public Test {
     public:
       /// Create and register test
       MinXXX(const std::string& s, const Gecode::IntSet& d,
              Gecode::IntConLevel icl)
         : Test("Arithmetic::Min::Bin::XXX::"+str(icl)+"::"+s,1,d) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::min(x[0],x[0]) == x[0];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::min(home, x[0], x[0], x[0], icl);
       }
     };

     /// %Test for binary maximum constraint
     class MaxXYZ : public Test {
     public:
       /// Create and register test
       MaxXYZ(const std::string& s, const Gecode::IntSet& d,
              Gecode::IntConLevel icl)
         : Test("Arithmetic::Max::Bin::XYZ::"+str(icl)+"::"+s,3,d) {
         contest = CTL_BOUNDS_Z;
       }
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::max(x[0],x[1]) == x[2];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::max(home, x[0], x[1], x[2], icl);
       }
     };

     /// %Test for binary maximum constraint with shared variables
     class MaxXXY : public Test {
     public:
       /// Create and register test
       MaxXXY(const std::string& s, const Gecode::IntSet& d,
              Gecode::IntConLevel icl)
         : Test("Arithmetic::Max::Bin::XXY::"+str(icl)+"::"+s,2,d) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::max(x[0],x[0]) == x[1];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::max(home, x[0], x[0], x[1], icl);
       }
     };

     /// %Test for binary maximum constraint with shared variables
     class MaxXYX : public Test {
     public:
       /// Create and register test
       MaxXYX(const std::string& s, const Gecode::IntSet& d,
              Gecode::IntConLevel icl)
         : Test("Arithmetic::Max::Bin::XYX::"+str(icl)+"::"+s,2,d) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::max(x[0],x[1]) == x[0];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::max(home, x[0], x[1], x[0], icl);
       }
     };

     /// %Test for binary maximum constraint with shared variables
     class MaxXYY : public Test {
     public:
       /// Create and register test
       MaxXYY(const std::string& s, const Gecode::IntSet& d,
              Gecode::IntConLevel icl)
         : Test("Arithmetic::Max::Bin::XYY::"+str(icl)+"::"+s,2,d) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::max(x[0],x[1]) == x[1];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::max(home, x[0], x[1], x[1], icl);
       }
     };

     /// %Test for binary maximum constraint with shared variables
     class MaxXXX : public Test {
     public:
       /// Create and register test
       MaxXXX(const std::string& s, const Gecode::IntSet& d,
              Gecode::IntConLevel icl)
         : Test("Arithmetic::Max::Bin::XXX::"+str(icl)+"::"+s,1,d) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::max(x[0],x[0]) == x[0];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::max(home, x[0], x[0], x[0], icl);
       }
     };

     /// %Test for n-ary minimmum constraint
     class MinNary : public Test {
     public:
       /// Create and register test
       MinNary(Gecode::IntConLevel icl)
         : Test("Arithmetic::Min::Nary::"+str(icl),4,-4,4,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::min(std::min(x[0],x[1]), x[2]) == x[3];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::IntVarArgs m(3);
         m[0]=x[0]; m[1]=x[1]; m[2]=x[2];
         Gecode::min(home, m, x[3], icl);
       }
     };

     /// %Test for n-ary minimmum constraint with shared variables
     class MinNaryShared : public Test {
     public:
       /// Create and register test
       MinNaryShared(Gecode::IntConLevel icl)
         : Test("Arithmetic::Min::Nary::Shared::"+str(icl),3,-4,4,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::min(std::min(x[0],x[1]), x[2]) == x[1];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::IntVarArgs m(3);
         m[0]=x[0]; m[1]=x[1]; m[2]=x[2];
         Gecode::min(home, m, x[1], icl);
       }
     };

     /// %Test for n-ary maximum constraint
     class MaxNary : public Test {
     public:
       /// Create and register test
       MaxNary(Gecode::IntConLevel icl)
         : Test("Arithmetic::Max::Nary::"+str(icl),4,-4,4,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::max(std::max(x[0],x[1]), x[2]) == x[3];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::IntVarArgs m(3);
         m[0]=x[0]; m[1]=x[1]; m[2]=x[2];
         Gecode::max(home, m, x[3], icl);
       }
     };

     /// %Test for n-ary maximum constraint with shared variables
     class MaxNaryShared : public Test {
     public:
       /// Create and register test
       MaxNaryShared(Gecode::IntConLevel icl)
         : Test("Arithmetic::Max::Nary::Shared::"+str(icl),3,-4,4,false,icl) {}
       /// %Test whether \a x is solution
       virtual bool solution(const Assignment& x) const {
         return std::max(std::max(x[0],x[1]), x[2]) == x[1];
       }
       /// Post constraint on \a x
       virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) {
         Gecode::IntVarArgs m(3);
         m[0]=x[0]; m[1]=x[1]; m[2]=x[2];
         Gecode::max(home, m, x[1], icl);
       }
     };


     /// Help class to create and register tests
     class Create {
     public:
       /// Perform creation and registration
       Create(void) {

         const int va[7] = {
           Gecode::Int::Limits::min, Gecode::Int::Limits::min+1,
           -1,0,1,
           Gecode::Int::Limits::max-1, Gecode::Int::Limits::max
         };
         const int vb[9] = {
           static_cast<int>(-sqrt(static_cast<double>
                                  (-Gecode::Int::Limits::min))),
           -4,-2,-1,0,1,2,4,
           static_cast<int>(sqrt(static_cast<double>
                                 (Gecode::Int::Limits::max)))
         };
         
         Gecode::IntSet a(va,7);
         Gecode::IntSet b(vb,9);
         Gecode::IntSet c(-8,8);
         Gecode::IntSet d(-70,70);

         (void) new DivMod("A",a);
         (void) new DivMod("B",b);
         (void) new DivMod("C",c);
         
         (void) new Div("A",a);
         (void) new Div("B",b);
         (void) new Div("C",c);
         
         (void) new Mod("A",a);
         (void) new Mod("B",b);
         (void) new Mod("C",c);


         for (IntConLevels icls; icls(); ++icls) 
           if (icls.icl() != Gecode::ICL_VAL) {
             (void) new MultXYZ("A",a,icls.icl());
             (void) new MultXYZ("B",b,icls.icl());
             (void) new MultXYZ("C",c,icls.icl());

             (void) new MultXXY("A",a,icls.icl());
             (void) new MultXXY("B",b,icls.icl());
             (void) new MultXXY("C",c,icls.icl());

             (void) new MultXYX("A",a,icls.icl());
             (void) new MultXYX("B",b,icls.icl());
             (void) new MultXYX("C",c,icls.icl());

             (void) new MultXYY("A",a,icls.icl());
             (void) new MultXYY("B",b,icls.icl());
             (void) new MultXYY("C",c,icls.icl());

             (void) new MultXXX("A",a,icls.icl());
             (void) new MultXXX("B",b,icls.icl());
             (void) new MultXXX("C",c,icls.icl());

             (void) new SqrXY("A",a,icls.icl());
             (void) new SqrXY("B",b,icls.icl());
             (void) new SqrXY("C",c,icls.icl());

             (void) new SqrXX("A",a,icls.icl());
             (void) new SqrXX("B",b,icls.icl());
             (void) new SqrXX("C",c,icls.icl());

             for (int n=0; n<=6; n++) {
               (void) new PowXY("A",n,a,icls.icl());
               (void) new PowXY("B",n,b,icls.icl());
               (void) new PowXY("C",n,c,icls.icl());
               (void) new PowXY("D",n,d,icls.icl());
  
               (void) new PowXX("A",n,a,icls.icl());
               (void) new PowXX("B",n,b,icls.icl());
               (void) new PowXX("C",n,c,icls.icl());
               (void) new PowXX("D",n,d,icls.icl());
             }

             for (int n=1; n<=6; n++) {
               (void) new NrootXY("A",n,a,icls.icl());
               (void) new NrootXY("B",n,b,icls.icl());
               (void) new NrootXY("C",n,c,icls.icl());
               (void) new NrootXY("D",n,d,icls.icl());
  
               (void) new NrootXX("A",n,a,icls.icl());
               (void) new NrootXX("B",n,b,icls.icl());
               (void) new NrootXX("C",n,c,icls.icl());
               (void) new NrootXX("D",n,d,icls.icl());
             }

             for (int n=30; n<=34; n++) {
               (void) new PowXY("C",n,c,icls.icl());
               (void) new PowXX("C",n,c,icls.icl());
               (void) new NrootXY("C",n,c,icls.icl());
               (void) new NrootXX("C",n,c,icls.icl());
             }

             (void) new SqrtXY("A",a,icls.icl());
             (void) new SqrtXY("B",b,icls.icl());
             (void) new SqrtXY("C",c,icls.icl());

             (void) new SqrtXX("A",a,icls.icl());
             (void) new SqrtXX("B",b,icls.icl());
             (void) new SqrtXX("C",c,icls.icl());

             (void) new AbsXY("A",a,icls.icl());
             (void) new AbsXY("B",b,icls.icl());
             (void) new AbsXY("C",c,icls.icl());

             (void) new AbsXX("A",a,icls.icl());
             (void) new AbsXX("B",b,icls.icl());
             (void) new AbsXX("C",c,icls.icl());

             (void) new MinXYZ("A",a,icls.icl());
             (void) new MinXYZ("B",b,icls.icl());
             (void) new MinXYZ("C",c,icls.icl());

             (void) new MinXXY("A",a,icls.icl());
             (void) new MinXXY("B",b,icls.icl());
             (void) new MinXXY("C",c,icls.icl());

             (void) new MinXYX("A",a,icls.icl());
             (void) new MinXYX("B",b,icls.icl());
             (void) new MinXYX("C",c,icls.icl());

             (void) new MinXYY("A",a,icls.icl());
             (void) new MinXYY("B",b,icls.icl());
             (void) new MinXYY("C",c,icls.icl());

             (void) new MinXXX("A",a,icls.icl());
             (void) new MinXXX("B",b,icls.icl());
             (void) new MinXXX("C",c,icls.icl());

             (void) new MaxXYZ("A",a,icls.icl());
             (void) new MaxXYZ("B",b,icls.icl());
             (void) new MaxXYZ("C",c,icls.icl());

             (void) new MaxXXY("A",a,icls.icl());
             (void) new MaxXXY("B",b,icls.icl());
             (void) new MaxXXY("C",c,icls.icl());

             (void) new MaxXYX("A",a,icls.icl());
             (void) new MaxXYX("B",b,icls.icl());
             (void) new MaxXYX("C",c,icls.icl());

             (void) new MaxXYY("A",a,icls.icl());
             (void) new MaxXYY("B",b,icls.icl());
             (void) new MaxXYY("C",c,icls.icl());

             (void) new MaxXXX("A",a,icls.icl());
             (void) new MaxXXX("B",b,icls.icl());
             (void) new MaxXXX("C",c,icls.icl());

             (void) new MinNary(icls.icl());
             (void) new MinNaryShared(icls.icl());
             (void) new MaxNary(icls.icl());
             (void) new MaxNaryShared(icls.icl());
           }
       }
     };

     Create c;
     //@}

   }
}}

// STATISTICS: test-int
