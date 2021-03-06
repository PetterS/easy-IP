/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Christian Schulte <schulte@gecode.org>
 *
 *  Copyright:
 *     Christian Schulte, 2009
 *
 *  Last modified:
 *     $Date: 2012-09-27 08:00:33 +0200 (to, 27 sep 2012) $ by $Author: tack $
 *     $Revision: 13118 $
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

#include <gecode/support.hh>

namespace Gecode { namespace Support {

  /*
   * Threads
   */
  
  Mutex* Thread::m(void) {
    static Mutex* m = new Mutex;
    return m;
  }
  
  Thread::Run* Thread::idle = NULL;

  void
  Thread::Run::exec(void) {
    while (true) {
      // Execute runnable
      {
        Runnable* e;
        m.acquire();
        e=r; r=NULL;
        m.release();
        assert(e != NULL);
        e->run();
        delete e;
      }
      // Put into idle stack
      Thread::m()->acquire();
      n=Thread::idle; Thread::idle=this;
      Thread::m()->release();
      // Wait for next runnable
      e.wait();
    }
  }

}}

// STATISTICS: support-any
