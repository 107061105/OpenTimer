#ifndef OT_TIMER_PFXT_HPP_
#define OT_TIMER_PFXT_HPP_

#include <ot/headerdef.hpp>

namespace ot {

// Forward declaration
class SfxtCache;
class Arc;

// ------------------------------------------------------------------------------------------------

// Class: PfxtNode
struct PfxtNode {
  
  Statisical_delay slack;
  size_t from;
  size_t to;

  const Arc* arc {nullptr};
  const PfxtNode* parent {nullptr};

  PfxtNode(Statisical_delay, size_t, size_t, const Arc*, const PfxtNode*);
};

// ------------------------------------------------------------------------------------------------

// Class: PfxtCache
class PfxtCache {

  friend class Timer;
  
  // Min-heap comparator
  struct PfxtNodeComparator {
    bool operator () (std::unique_ptr<PfxtNode>& a, std::unique_ptr<PfxtNode>& b) const {
      return a->slack.mean() > b->slack.mean(); // TODO XD
    }
  };

  public:
    
    inline size_t num_nodes() const;

    PfxtCache(const SfxtCache&);

    PfxtCache(const PfxtCache&) = delete;
    PfxtCache(PfxtCache&&);

    PfxtCache& operator = (const PfxtCache&) = delete;
    PfxtCache& operator = (PfxtCache&&) = delete;

  private:

    const SfxtCache& _sfxt;

    PfxtNodeComparator _comp;

    std::vector<std::unique_ptr<PfxtNode>> _paths;
    std::vector<std::unique_ptr<PfxtNode>> _nodes;

    void _push(Statisical_delay, size_t, size_t, const Arc*, const PfxtNode*);

    PfxtNode* _pop();
    PfxtNode* _top() const;
};

// Function: num_nodes
inline size_t PfxtCache::num_nodes() const {
  return _nodes.size();
}

};  // end of namespace ot. ----------------------------------------------------------------------

#endif




