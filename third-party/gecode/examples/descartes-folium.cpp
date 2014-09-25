/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Vincent Barichard <Vincent.Barichard@univ-angers.fr>
 *
 *  Copyright:
 *     Vincent Barichard, 2012
 *
 *  Last modified:
 *     $Date: 2013-07-08 14:22:40 +0200 (må, 08 jul 2013) $ by $Author: schulte $
 *     $Revision: 13820 $
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

#include <gecode/driver.hh>

#include <gecode/minimodel.hh>
#include <gecode/float.hh>

using namespace Gecode;

/**
 * \brief %Example: Folium of Descartes
 *
 * The folium of Descartes is a curve defined by the equation:
 * \f[
 * x^3 + y^3 - 3axy = 0
 * \f]
 * 
 * A technique to solve it, is to write \f$y=px\f$ and solve for
 * \f$x\f$ and \f$y\f$ in terms of \f$p\f$.  By setting
 * \f$a=1\f$, it yields to the paramatric equation:
 * 
 * \f[
 * x^3 + y^3 - 3xy = 0
 * \f]
 * \f[
 * x=\frac{3p}{1+p^3},\quad y=\frac{3p^2}{1+p^3}
 * \f]
 * 
 * The parameter \f$p\f$ is related to the position on the curve
 * and is constrained to get different solutions for \f$x\f$ and
 * \f$y\f$. To get reasonable interval starting sizes, \f$p\f$
 * and \f$y\f$ are restricted to \f$[-20;20]\f$ and \f$x\f$ is
 * restricted to \f$[-1;2]\f$.
 *
 * \ingroup Example
 */
class DescartesFolium : public Script {
protected:
  /// The numbers
  FloatVarArray f;
  /// Minimum distance between two solutions
  double step;
public:
  /// Actual model
  DescartesFolium(const Options&) 
    : f(*this,3,-20,20), step(0.1) {
    // Post equation
    FloatVar p = f[0];
    FloatVar x = f[1];
    FloatVar y = f[2];
    rel(*this, 3*p/(1+pow(p,3)) == x);
    rel(*this, 3*sqr(p)/(1+pow(p,3)) == y);
    rel(*this, pow(x,3) + pow(y,3)  == 3 * x * y);
    rel(*this, x == FloatVal(-1,2));

    branch(*this,p,FLOAT_VAL_SPLIT_MIN());
  }
  /// Constructor for cloning \a p
  DescartesFolium(bool share, DescartesFolium& p) 
    : Script(share,p), step(p.step) {
    f.update(*this,share,p.f);
  }
  /// Copy during cloning
  virtual Space* copy(bool share) { 
    return new DescartesFolium(share,*this); 
  }
  /// Add constraint to current model to get next solution (not too close)
  virtual void constrain(const Space& _b) {
    const DescartesFolium& b = static_cast<const DescartesFolium&>(_b);
    rel(*this, f[0] >= (b.f[0].max()+step));
  }
  /// Print solution coordinates
  virtual void print(std::ostream& os) const {
    os << "XY " << f[1].med() << " " << f[2].med()
       << std::endl;
  }

};

/** \brief Main-function
 *  \relates DescartesFolium
 */
int main(int argc, char* argv[]) {
  Options opt("DescartesFolium");
  opt.parse(argc,argv);
  opt.solutions(0);
  Script::run<DescartesFolium,BAB,Options>(opt);
  return 0;
}

// STATISTICS: example-any
