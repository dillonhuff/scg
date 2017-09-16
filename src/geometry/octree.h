//
//  Octree.h
//
//  Created by Eduardo Poyart on 6/4/12.
//

/*
Copyright (c) 2012, Eduardo Poyart.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the
distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef GCA_OCTREE_H
#define GCA_OCTREE_H

#include <assert.h>

#define COMPUTE_SIDE(i, bit, p, mid, newMin, newMax) \
if (p >= mid)         \
{                     \
    i |= bit;         \
    newMin = mid;     \
}                     \
else                  \
{                     \
    newMax = mid;     \
}


namespace gca {

  template <class N>
  class octree
  {
  protected:
    struct Point
    {
      double x; 
      double y; 
      double z;
      Point(const Point& p2): x(p2.x), y(p2.y), z(p2.z) {}
      Point& operator=(const Point& p2) { x = p2.x; y = p2.y; z = p2.z; return *this;}
      Point(double in_x, double in_y, double in_z): x(in_x), y(in_y), z(in_z) {}
      Point(const double p2[3]): x(p2[0]), y(p2[1]), z(p2[2]) {}
      operator double*() { return &x; }
      operator const double*() const { return &x; }
      Point operator+(const Point& p2) const { return Point(x+p2.x, y+p2.y, z+p2.z); }
      Point operator-(const Point& p2) const { return Point(x-p2.x, y-p2.y, z-p2.z); }
      Point operator*(double f) const { return Point(x*f, y*f, z*f); }
      bool operator< (const Point& p2) const { return x <  p2.x && y <  p2.y && z <  p2.z; }
      bool operator>=(const Point& p2) const { return x >= p2.x && y >= p2.y && z >= p2.z; }
    };
    
    struct octree_node
    {
      N _nodeData;
      octree_node* _children[8];
      octree_node()
      {
	for (int i = 0; i < 8; i++)
	  _children[i] = 0;
      }
      virtual ~octree_node()
      {
	for (int i = 0; i < 8; i++)
	  if (_children[i])
	    delete _children[i];
      }
    };
    
    Point _min;
    Point _max;
    Point _cellSize;
    octree_node* _root;
    
  public:
    octree(double min[3], double max[3], double cellSize[3]): _min(min), _max(max), _cellSize(cellSize), _root(0) {}
    virtual ~octree() { delete _root; }
    
    class Callback
    {
    public:
      // Return value: true = continue; false = abort.
      virtual bool operator()(const double min[3], const double max[3], N& nodeData) = 0;
    };
    
    N& getCell(const double pos[3], Callback* callback = NULL)
    {
      Point ppos(pos);
      assert(ppos >= _min && ppos < _max);
      Point currMin(_min);
      Point currMax(_max);
      Point delta = _max - _min;
      if (!_root)
	_root = new octree_node();
      octree_node* currnode = _root;
      while (delta >= _cellSize)
        {
	  bool shouldContinue = true;
	  if (callback)
	    shouldContinue = callback->operator()(currMin, currMax, currnode->_nodeData);
	  if (!shouldContinue)
	    break;
	  Point mid = (delta * 0.5f) + currMin;
	  Point newMin(currMin);
	  Point newMax(currMax);
	  int index = 0;
	  COMPUTE_SIDE(index, 1, ppos.x, mid.x, newMin.x, newMax.x)
            COMPUTE_SIDE(index, 2, ppos.y, mid.y, newMin.y, newMax.y)
            COMPUTE_SIDE(index, 4, ppos.z, mid.z, newMin.z, newMax.z)
            if (!(currnode->_children[index]))
	      currnode->_children[index] = new octree_node();
	  currnode = currnode->_children[index];
	  currMin = newMin;
	  currMax = newMax;
	  delta = currMax - currMin;
        }
      return currnode->_nodeData;
    }
    
    void traverse(Callback* callback)
    {
      assert(callback);
      traverseRecursive(callback, _min, _max, _root);
    }
    
    void clear()
    {
      delete _root;
      _root = NULL;
    }
    
    class Iterator {
    public:
      Iterator getChild(int i)
      {
	return Iterator(_currnode->_children[i]);
      }
      N* getData()
      {
	if (_currnode)
	  return &_currnode->_nodeData;
	else return NULL;
      }
    protected:
      octree_node* _currnode;
      Iterator(octree_node* node): _currnode(node) {}
      friend class octree;
    };
    
    Iterator getIterator()
    {
      return Iterator(_root);
    }
    
  protected:
    void traverseRecursive(Callback* callback, const Point& currMin, const Point& currMax, octree_node* currnode)
    {
      if (!currnode)
	return;
      bool shouldContinue = callback->operator()(currMin, currMax, currnode->_nodeData);
      if (!shouldContinue)
	return;
      Point delta = currMax - currMin;
      Point mid = (delta * 0.5f) + currMin;
      traverseRecursive(callback, currMin, mid, currnode->_children[0]);
      traverseRecursive(callback, Point(mid.x, currMin.y, currMin.z), 
			Point(currMax.x, mid.y, mid.z), currnode->_children[1]);
      traverseRecursive(callback, Point(currMin.x, mid.y, currMin.z), 
			Point(mid.x, currMax.y, mid.z), currnode->_children[2]);
      traverseRecursive(callback, Point(mid.x, mid.y, currMin.z), 
			Point(currMax.x, currMax.y, mid.z), currnode->_children[3]);
      traverseRecursive(callback, Point(currMin.x, currMin.y, mid.z), 
			Point(mid.x, mid.y, currMax.z), currnode->_children[4]);
      traverseRecursive(callback, Point(mid.x, currMin.y, mid.z), 
			Point(currMax.x, mid.y, currMax.z), currnode->_children[5]);
      traverseRecursive(callback, Point(currMin.x, mid.y, mid.z), 
			Point(mid.x, currMax.y, currMax.z), currnode->_children[6]);
      traverseRecursive(callback, mid, currMax, currnode->_children[7]);
    }    
  };

}
#endif
