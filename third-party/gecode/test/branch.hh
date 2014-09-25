/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Mikael Lagerkvist <lagerkvist@gecode.org>
 *     Christian Schulte <schulte@gecode.org>
 *
 *  Contributing authors:
 *     Vincent Barichard <Vincent.Barichard@univ-angers.fr>
 *
 *  Copyright:
 *     Mikael Lagerkvist, 2006
 *     Christian Schulte, 2009
 *     Vincent Barichard, 2012
 *
 *  Last modified:
 *     $Date: 2012-04-05 12:00:11 +0200 (to, 05 apr 2012) $ by $Author: vbarichard $
 *     $Revision: 12703 $
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

#ifndef __GECODE_TEST_BRANCH_HH__
#define __GECODE_TEST_BRANCH_HH__

#include <gecode/kernel.hh>

#include <gecode/int.hh>

#ifdef GECODE_HAS_SET_VARS
#include <gecode/set.hh>
#endif

#ifdef GECODE_HAS_FLOAT_VARS
#include <gecode/float.hh>
#endif

#include "test/test.hh"

namespace Test {

  /// Tests for branchings
  namespace Branch {

    /**
     * \brief %Base class for tests for branching on integer variables
     *
     */
    class IntTest : public Base {
    protected:
      /// Number of variables
      int arity;
      /// Domain of variables
      Gecode::IntSet dom;
    public:
      /// Construct and register test
      IntTest(const std::string& s, int a, const Gecode::IntSet& d);
      /// Perform test
      virtual bool run(void);
      /// Post propagators on variables \a x
      virtual void post(Gecode::Space& home, Gecode::IntVarArray& x) = 0;
    };

    /**
     * \brief %Base class for tests for branching on Boolean variables
     *
     */
    class BoolTest : public Base {
    protected:
      /// Number of variables
      int arity;
    public:
      /// Construct and register test
      BoolTest(const std::string& s, int a);
      /// Perform test
      virtual bool run(void);
      /// Post propagators on variables \a x
      virtual void post(Gecode::Space& home, Gecode::BoolVarArray& x) = 0;
    };

#ifdef GECODE_HAS_SET_VARS
    /**
     * \brief %Base class for tests for branching on set variables
     *
     */
    class SetTest : public Base {
    protected:
      /// Number of variables
      int arity;
      /// Domain of variables
      Gecode::IntSet dom;
    public:
      /// Construct and register test
      SetTest(const std::string& s, int a, const Gecode::IntSet& d);
      /// Perform test
      virtual bool run(void);
      /// Post propagators on variables \a x
      virtual void post(Gecode::Space& home, Gecode::SetVarArray& x) = 0;
    };
#endif

#ifdef GECODE_HAS_FLOAT_VARS
    /**
     * \brief %Base class for tests for branching on float variables
     *
     */
    class FloatTest : public Base {
    protected:
      /// Number of variables
      int arity;
      /// Domain of variables
      Gecode::FloatVal dom;
      /// Maximum number of solutions searched during solving
      int nbSols;
    public:
      /// Construct and register test
      FloatTest(const std::string& s, int a, const Gecode::FloatVal& d, int nbs);
      /// Perform test
      virtual bool run(void);
      /// Post propagators on variables \a x
      virtual void post(Gecode::Space& home, Gecode::FloatVarArray& x) = 0;
    };
#endif

  }

}

#endif

// STATISTICS: test-branch
