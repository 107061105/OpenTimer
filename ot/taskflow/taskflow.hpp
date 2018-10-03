// MIT License
// 
// Copyright (c) 2018 Tsung-Wei Huang, Chun-Xun Lin, Guannan Guo, and Martin Wong
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <iostream>
#include <mutex>
#include <deque>
#include <vector>
#include <algorithm>
#include <thread>
#include <future>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <list>
#include <forward_list>
#include <numeric>
#include <iomanip>
#include <cassert>
#include <optional>

#include "threadpool/threadpool.hpp"

// ============================================================================

// Clang mis-interprets variant's get as a non-friend of variant and cannot
// get compiled correctly. We use the patch: 
// https://gcc.gnu.org/viewcvs/gcc?view=revision&revision=258854
// to get rid of this.
#if defined(__clang__)
  #include "patch/clang_variant.hpp"
#else
  #include <variant>
#endif

// ============================================================================

/*// Class: ObjectPool
template <typename T>
class ObjectPool {
  
  struct Deleter {
    Deleter(std::vector<T*>&);
    void operator()(T*);
    std::vector<T*>& recycle;
  };
  
  public:
  
  using HandleType = std::unique_ptr<T, Deleter>;
    
    template <typename... ArgsT>
    auto get(ArgsT&&...);

  private:
  
    std::forward_list<T> _pool;
    std::vector<T*> _recycle; 
};

// Constructor
template <typename T>
ObjectPool<T>::Deleter::Deleter(std::vector<T*>& r) : recycle(r) {
}

// Operator
template <typename T>
void ObjectPool<T>::Deleter::operator()(T* item) {
  if(item != nullptr) {
    item->~T();
    recycle.push_back(item);
  }
}

// Constructor
template <typename T>
template <typename... ArgsT>
auto ObjectPool<T>::get(ArgsT&&... args) {
  // Pool is full
  if(_recycle.empty()) {
    T& item = _pool.emplace_front(std::forward<ArgsT>(args)...);
    return HandleType(&item, Deleter(_recycle));
  }
  // Get item from the recycle box
  else {
    auto item = _recycle.back(); 
    _recycle.pop_back();
    new (item) T(std::forward<ArgsT>(args)...);
    return HandleType(item, Deleter(_recycle));
  }
}*/


// Namespace of taskflow. -----------------------------------------------------
namespace tf {

// Procedure: throw_re
template <typename... ArgsT>
inline void throw_re(const char* fname, const size_t line, ArgsT&&... args) {
  std::ostringstream oss;
  oss << '[' << fname << ':' << line << "] ";
  (oss << ... << std::forward<ArgsT>(args));
  throw std::runtime_error(oss.str());
}

#define TF_THROW(...) throw_re(__FILE__, __LINE__, __VA_ARGS__);

//-----------------------------------------------------------------------------
// Traits
//-----------------------------------------------------------------------------

// Macro to check whether a class has a member function
#define define_has_member(member_name)                                     \
template <typename T>                                                      \
class has_member_##member_name                                             \
{                                                                          \
  typedef char yes_type;                                                   \
  typedef long no_type;                                                    \
  template <typename U> static yes_type test(decltype(&U::member_name));   \
  template <typename U> static no_type  test(...);                         \
  public:                                                                  \
    static constexpr bool value = sizeof(test<T>(0)) == sizeof(yes_type);  \
}

#define has_member(class_, member_name)  has_member_##member_name<class_>::value

// Struct: dependent_false
template <typename... T>
struct dependent_false { 
  static constexpr bool value = false; 
};

template <typename... T>
constexpr auto dependent_false_v = dependent_false<T...>::value;

// Struct: is_iterator
template <typename T, typename = void>
struct is_iterator {
  static constexpr bool value = false;
};

template <typename T>
struct is_iterator<
  T, 
  std::enable_if_t<!std::is_same_v<typename std::iterator_traits<T>::value_type, void>>
> {
  static constexpr bool value = true;
};

template <typename T>
inline constexpr bool is_iterator_v = is_iterator<T>::value;

// Struct: is_iterable
template <typename T, typename = void>
struct is_iterable : std::false_type {
};

template <typename T>
struct is_iterable<T, std::void_t<decltype(std::declval<T>().begin()),
                                  decltype(std::declval<T>().end())>>
  : std::true_type {
};

template <typename T>
inline constexpr bool is_iterable_v = is_iterable<T>::value;

//-----------------------------------------------------------------------------
// Taskflow definition
//-----------------------------------------------------------------------------



// Forward declaration
template <template<typename...> typename Func>
class BasicNode;

template <typename Node>
class BasicTopology;

template <typename Node>
class BasicTask;

template <typename Node>
class BasicFlowBuilder;

template <typename Node>
class BasicSubflowBuilder;

template <template <typename...> typename F, typename Threadpool>
class BasicTaskflow;

// ----------------------------------------------------------------------------

// Class: BasicNode
template <template<typename...> typename F>
class BasicNode {

  template <typename U> friend class BasicTask;
  template <typename S> friend class BasicTopology;
  template <template <typename...> class X, typename Y> friend class BasicTaskflow;

  using Work     = F<void()>;
  using Subwork  = F<void(BasicSubflowBuilder<BasicNode>&)>;
  using Topology = BasicTopology<BasicNode>;

  public:

    BasicNode();

    template <typename C>
    BasicNode(C&&);

    const std::string& name() const;
    
    void precede(BasicNode&);

    size_t num_successors() const;
    size_t num_dependents() const;

    std::string dump() const;

  private:
    
    std::string _name;
    std::variant<Work, Subwork> _work;
    std::vector<BasicNode*> _successors;
    std::atomic<int> _dependents;
    std::forward_list<BasicNode> _children;
    Topology* _topology;

    void _dump(std::ostream&) const;
};

// Constructor
template <template<typename...> typename F>
BasicNode<F>::BasicNode() {
  _dependents.store(0, std::memory_order_relaxed);
  _topology = nullptr;
}

// Constructor
template <template<typename...> typename F>
template <typename C>
BasicNode<F>::BasicNode(C&& c) : _work {std::forward<C>(c)} {
  _dependents.store(0, std::memory_order_relaxed);
  _topology = nullptr;
}

// Procedure:
template <template<typename...> typename F>
void BasicNode<F>::precede(BasicNode& v) {
  _successors.push_back(&v);
  //++v._dependents;
  v._dependents.fetch_add(1, std::memory_order_relaxed);
}

// Function: num_successors
template <template<typename...> typename F>
size_t BasicNode<F>::num_successors() const {
  return _successors.size();
}

// Function: dependents
template <template<typename...> typename F>
size_t BasicNode<F>::num_dependents() const {
  return _dependents.load(std::memory_order_relaxed);
}

// Function: name
template <template<typename...> typename F>
const std::string& BasicNode<F>::name() const {
  return _name;
}

// Function: dump
template <template<typename...> typename F>
std::string BasicNode<F>::dump() const {
  std::ostringstream os;  
  _dump(os);
  return os.str();
}

// Function: _dump
template <template<typename...> typename F>
void BasicNode<F>::_dump(std::ostream& os) const {
  
  if(_name.empty()) os << '\"' << this << '\"';
  else os << std::quoted(_name);
  os << ";\n";

  for(const auto s : _successors) {

    if(_name.empty()) os << '\"' << this << '\"';
    else os << std::quoted(_name);

    os << " -> ";
    
    if(s->name().empty()) os << '\"' << s << '\"';
    else os << std::quoted(s->name());

    os << ";\n";
  }
  
  if(!_children.empty()) {

    os << "subgraph cluster_";
    if(_name.empty()) os << this;
    else os << _name;
    os << " {\n";

    os << "label = \"Subflow_";
    if(_name.empty()) os << this;
    else os << _name;

    os << "\";\n" << "color=blue\n";

    for(const auto& n : _children) {
      n._dump(os);
    }
    os << "}\n";
  }
}

// ----------------------------------------------------------------------------
  
// class: BasicTopology
template <typename Node>
class BasicTopology {
  
  template <template <typename...> class F, typename P> friend class BasicTaskflow;

  public:

    BasicTopology(std::forward_list<Node>&&);

    std::string dump() const;

  private:

    std::forward_list<Node> _nodes;
    std::shared_future<void> _future;

    Node _source;
    Node _target;

    void _dump(std::ostream&) const;
};

// Constructor
template <typename Node>
BasicTopology<Node>::BasicTopology(std::forward_list<Node>&& t) : 
  _nodes(std::move(t)) {

  _source._topology = this;
  _target._topology = this;
  
  std::promise<void> promise;

  _future = promise.get_future().share();

  _target._work = [p=MoC{std::move(promise)}] () mutable { 
    p.get().set_value(); 
  };

  // ensure the topology is connected
  _source.precede(_target);

  // Build the super source and super target.
  for(auto& node : _nodes) {

    node._topology = this;

    if(node.num_dependents() == 0) {
      _source.precede(node);
    }

    if(node.num_successors() == 0) {
      node.precede(_target);
    }
  }
}

// Procedure: _dump
template <typename Node>
void BasicTopology<Node>::_dump(std::ostream& os) const {

  assert(_source._children.empty());
  assert(_target._children.empty());
  
  os << "digraph Topology {\n"
     << _source.dump() 
     << _target.dump();

  for(const auto& node : _nodes) {
    os << node.dump();
  }

  os << "}\n";
}
  
// Function: dump
template <typename Node>
std::string BasicTopology<Node>::dump() const { 
  std::ostringstream os;
  _dump(os);
  return os.str();
}

// ----------------------------------------------------------------------------

// Class: BasicTask
template <typename Node>
class BasicTask {

  template <typename U> friend class BasicFlowBuilder;
  template <template <typename...> class F, typename P> friend class BasicTaskflow;

  public:
    
    BasicTask() = default;
    BasicTask(Node&);
    BasicTask(const BasicTask&);
    BasicTask(BasicTask&&);

    BasicTask& operator = (const BasicTask&);

    const std::string& name() const;

    size_t num_successors() const;
    size_t num_dependents() const;

    BasicTask& name(const std::string&);
    BasicTask& precede(BasicTask);
    BasicTask& broadcast(std::vector<BasicTask>&);
    BasicTask& broadcast(std::initializer_list<BasicTask>);
    BasicTask& gather(std::vector<BasicTask>&);
    BasicTask& gather(std::initializer_list<BasicTask>);

    template <typename C>
    BasicTask& work(C&&);
  
    template <typename... Bs>
    BasicTask& broadcast(Bs&&...);

    template <typename... Bs>
    BasicTask& gather(Bs&&...);

  private:

    Node* _node {nullptr};

    template<typename S>
    void _broadcast(S&);

    template<typename S>
    void _gather(S&);
};

// Constructor
template <typename Node>
BasicTask<Node>::BasicTask(Node& t) : _node {&t} {
}

// Constructor
template <typename Node>
BasicTask<Node>::BasicTask(const BasicTask& rhs) : _node {rhs._node} {
}

// Function: precede
template <typename Node>
BasicTask<Node>& BasicTask<Node>::precede(BasicTask tgt) {
  _node->precede(*(tgt._node));
  return *this;
}

// Function: broadcast
template <typename Node>
template <typename... Bs>
BasicTask<Node>& BasicTask<Node>::broadcast(Bs&&... tgts) {
  (_node->precede(*(tgts._node)), ...);
  return *this;
}

// Procedure: _broadcast
template <typename Node>
template <typename S>
void BasicTask<Node>::_broadcast(S& tgts) {
  for(auto& to : tgts) {
    _node->precede(*(to._node));
  }
}
      
// Function: broadcast
template <typename Node>
BasicTask<Node>& BasicTask<Node>::broadcast(std::vector<BasicTask>& tgts) {
  _broadcast(tgts);
  return *this;
}

// Function: broadcast
template <typename Node>
BasicTask<Node>& BasicTask<Node>::broadcast(
  std::initializer_list<BasicTask> tgts
) {
  _broadcast(tgts);
  return *this;
}

// Function: gather
template <typename Node>
template <typename... Bs>
BasicTask<Node>& BasicTask<Node>::gather(Bs&&... tgts) {
  (tgts.precede(*this), ...);
  return *this;
}

// Procedure: _gather
template <typename Node>
template <typename S>
void BasicTask<Node>::_gather(S& tgts) {
  for(auto& from : tgts) {
    from._node->precede(*_node);
  }
}

// Function: gather
template <typename Node>
BasicTask<Node>& BasicTask<Node>::gather(std::vector<BasicTask>& tgts) {
  _gather(tgts);
  return *this;
}

// Function: gather
template <typename Node>
BasicTask<Node>& BasicTask<Node>::gather(
  std::initializer_list<BasicTask> tgts
) {
  _gather(tgts);
  return *this;
}

// Operator =
template <typename Node>
BasicTask<Node>& BasicTask<Node>::operator = (const BasicTask& rhs) {
  _node = rhs._node;
  return *this;
}

// Constructor
template <typename Node>
BasicTask<Node>::BasicTask(BasicTask&& rhs) : _node{rhs._node} { 
  rhs._node = nullptr; 
}

// Function: work
template <typename Node>
template <typename C>
BasicTask<Node>& BasicTask<Node>::work(C&& c) {
  _node->_work = std::forward<C>(c);
  return *this;
}

// Function: name
template <typename Node>
BasicTask<Node>& BasicTask<Node>::name(const std::string& name) {
  _node->_name = name;
  return *this;
}

// Function: name
template <typename Node>
const std::string& BasicTask<Node>::name() const {
  return _node->_name;
}

// Function: num_dependents
template <typename Node>
size_t BasicTask<Node>::num_dependents() const {
  return _node->_dependents.load(std::memory_order_relaxed);
}

// Function: num_successors
template <typename Node>
size_t BasicTask<Node>::num_successors() const {
  return _node->_successors.size();
}

// ----------------------------------------------------------------------------

// Class: BasicFlowBuilder
template <typename Node>
class BasicFlowBuilder {

  using Task = BasicTask<Node>;

  public:

    BasicFlowBuilder(std::forward_list<Node>&, size_t);

    template <typename C>
    auto emplace(C&&);

    template <typename... C, std::enable_if_t<(sizeof...(C)>1), void>* = nullptr>
    auto emplace(C&&...);

    template <typename C>
    auto silent_emplace(C&&);

    template <typename... C, std::enable_if_t<(sizeof...(C)>1), void>* = nullptr>
    auto silent_emplace(C&&...);

    template <typename I, typename C>
    auto parallel_for(I, I, C&&, size_t = 0);

    template <typename T, typename C, std::enable_if_t<is_iterable_v<T>, void>* = nullptr>
    auto parallel_for(T&, C&&, size_t = 0);

    template <typename I, typename T, typename B>
    auto reduce(I, I, T&, B&&);

    template <typename I, typename T>
    auto reduce_min(I, I, T&);
    
    template <typename I, typename T>
    auto reduce_max(I, I, T&);

    template <typename I, typename T, typename B, typename U>
    auto transform_reduce(I, I, T&, B&&, U&&);

    template <typename I, typename T, typename B, typename P, typename U>
    auto transform_reduce(I, I, T&, B&&, P&&, U&&);
    
    auto placeholder();
    
    void precede(Task, Task);
    void linearize(std::vector<Task>&);
    void linearize(std::initializer_list<Task>);
    void broadcast(Task, std::vector<Task>&);
    void broadcast(Task, std::initializer_list<Task>);
    void gather(std::vector<Task>&, Task);
    void gather(std::initializer_list<Task>, Task);  

    size_t num_nodes() const;

    bool empty() const;

  protected:

    std::forward_list<Node>& _nodes;
    size_t _num_workers;

    template <typename L>
    void _linearize(L&);
};

template <typename Node>    
BasicFlowBuilder<Node>::BasicFlowBuilder(
  std::forward_list<Node>& nodes, size_t num_workers
) : 
  _nodes       {nodes}, 
  _num_workers {num_workers} {
}    

// Procedure: num_nodes
template <typename Node>
size_t BasicFlowBuilder<Node>::num_nodes() const {
  return std::distance(_nodes.begin(), _nodes.end());
}

// Function: empty
template <typename Node>
bool BasicFlowBuilder<Node>::empty() const {
  return _nodes.empty();
}

// Procedure: precede
template <typename Node>
void BasicFlowBuilder<Node>::precede(Task from, Task to) {
  from._node->precede(*(to._node));
}

// Procedure: broadcast
template <typename Node>
void BasicFlowBuilder<Node>::broadcast(
  Task from, std::vector<Task>& keys
) {
  from.broadcast(keys);
}

// Procedure: broadcast
template <typename Node>
void BasicFlowBuilder<Node>::broadcast(
  Task from, std::initializer_list<Task> keys
) {
  from.broadcast(keys);
}

// Function: gather
template <typename Node>
void BasicFlowBuilder<Node>::gather(
  std::vector<Task>& keys, Task to
) {
  to.gather(keys);
}

// Function: gather
template <typename Node>
void BasicFlowBuilder<Node>::gather(
  std::initializer_list<Task> keys, Task to
) {
  to.gather(keys);
}

// Function: placeholder
template <typename Node>
auto BasicFlowBuilder<Node>::placeholder() {
  auto& node = _nodes.emplace_front();
  return Task(node);
}

// Function: emplace
template <typename Node>
template <typename C>
auto BasicFlowBuilder<Node>::emplace(C&& c) {
    
  // subflow task
  if constexpr(std::is_invocable_v<C, BasicSubflowBuilder<Node>&>) {

    using R = std::invoke_result_t<C, BasicSubflowBuilder<Node>&>;
    std::promise<R> p;
    auto fu = p.get_future();
  
    if constexpr(std::is_same_v<void, R>) {
      auto& node = _nodes.emplace_front([p=MoC(std::move(p)), c=std::forward<C>(c)]
      (BasicSubflowBuilder<Node>& fb) mutable {
        if(fb._nodes.empty()) {
          c(fb);
          if(fb.detached()) {
            p.get().set_value();
          }
        }
        else {
          p.get().set_value();
        }
      });
      return std::make_pair(Task(node), std::move(fu));
    }
    else {
      auto& node = _nodes.emplace_front(
      [p=MoC(std::move(p)), c=std::forward<C>(c), r=std::optional<R>()]
      (BasicSubflowBuilder<Node>& fb) mutable {
        if(fb._nodes.empty()) {
          r.emplace(c(fb));
          if(fb.detached()) {
            p.get().set_value(std::move(*r)); 
          }
        }
        else {
          assert(r);
          p.get().set_value(std::move(*r));
        }
      });
      return std::make_pair(Task(node), std::move(fu));
    }
  }
  // regular task
  else if constexpr(std::is_invocable_v<C>) {

    using R = std::invoke_result_t<C>;
    std::promise<R> p;
    auto fu = p.get_future();

    if constexpr(std::is_same_v<void, R>) {
      auto& node = _nodes.emplace_front(
        [p=MoC(std::move(p)), c=std::forward<C>(c)]() mutable {
          c(); 
          p.get().set_value();
        }
      );
      return std::make_pair(Task(node), std::move(fu));
    }
    else {
      auto& node = _nodes.emplace_front(
        [p=MoC(std::move(p)), c=std::forward<C>(c)]() mutable {
          p.get().set_value(c());
        }
      );
      return std::make_pair(Task(node), std::move(fu));
    }
  }
  else {
    static_assert(dependent_false_v<C>, "invalid task work type");
  }
}

// Function: emplace
template <typename Node>
template <typename... C, std::enable_if_t<(sizeof...(C)>1), void>*>
auto BasicFlowBuilder<Node>::emplace(C&&... cs) {
  return std::make_tuple(emplace(std::forward<C>(cs))...);
}

// Function: silent_emplace
template <typename Node>
template <typename C>
auto BasicFlowBuilder<Node>::silent_emplace(C&& c) {
  // subflow task
  if constexpr(std::is_invocable_v<C, BasicSubflowBuilder<Node>&>) {
    auto& n = _nodes.emplace_front(
    [c=std::forward<C>(c)] (BasicSubflowBuilder<Node>& fb) {
      // first time execution
      if(fb._nodes.empty()) {
        c(fb);
      }
    });
    return Task(n);
  }
  // regular task
  else if constexpr(std::is_invocable_v<C>) {
    auto& n = _nodes.emplace_front(std::forward<C>(c));
    return Task(n);
  }
  else {
    static_assert(dependent_false_v<C>, "invalid task work type");
  }
}

// Function: silent_emplace
template <typename Node>
template <typename... C, std::enable_if_t<(sizeof...(C)>1), void>*>
auto BasicFlowBuilder<Node>::silent_emplace(C&&... cs) {
  return std::make_tuple(silent_emplace(std::forward<C>(cs))...);
}


// Function: parallel_for    
template <typename Node>
template <typename I, typename C>
auto BasicFlowBuilder<Node>::parallel_for(I beg, I end, C&& c, size_t g) {

  using category = typename std::iterator_traits<I>::iterator_category;
  
  if(g == 0) {
    auto d = std::distance(beg, end);
    auto w = std::max(size_t{1}, _num_workers);
    g = (d + w - 1) / w;
  }

  auto source = placeholder();
  auto target = placeholder();
  
  while(beg != end) {

    auto e = beg;
    
    // Case 1: random access iterator
    if constexpr(std::is_same_v<category, std::random_access_iterator_tag>) {
      size_t r = std::distance(beg, end);
      std::advance(e, std::min(r, g));
    }
    // Case 2: non-random access iterator
    else {
      for(size_t i=0; i<g && e != end; ++e, ++i);
    }
      
    // Create a task
    auto task = silent_emplace([beg, e, c] () mutable {
      std::for_each(beg, e, c);
    });
    source.precede(task);
    task.precede(target);

    // adjust the pointer
    beg = e;
  }

  return std::make_pair(source, target); 
}

// Function: parallel_for
template <typename Node>
template <typename T, typename C, std::enable_if_t<is_iterable_v<T>, void>*>
auto BasicFlowBuilder<Node>::parallel_for(T& t, C&& c, size_t group) {
  return parallel_for(t.begin(), t.end(), std::forward<C>(c), group);
}

// Function: reduce_min
// Find the minimum element over a range of items.
template <typename Node>
template <typename I, typename T>
auto BasicFlowBuilder<Node>::reduce_min(I beg, I end, T& result) {
  return reduce(beg, end, result, [] (const auto& l, const auto& r) {
    return std::min(l, r);
  });
}

// Function: reduce_max
// Find the maximum element over a range of items.
template <typename Node>
template <typename I, typename T>
auto BasicFlowBuilder<Node>::reduce_max(I beg, I end, T& result) {
  return reduce(beg, end, result, [] (const auto& l, const auto& r) {
    return std::max(l, r);
  });
}

// Function: transform_reduce    
template <typename Node>
template <typename I, typename T, typename B, typename U>
auto BasicFlowBuilder<Node>::transform_reduce(
  I beg, I end, T& result, B&& bop, U&& uop
) {

  using category = typename std::iterator_traits<I>::iterator_category;
  
  // Even partition
  size_t d = std::distance(beg, end);
  size_t w = std::max(size_t{1}, _num_workers);
  size_t g = std::max((d + w - 1) / w, size_t{2});

  auto source = placeholder();
  auto target = placeholder();

  std::vector<std::future<T>> futures;

  while(beg != end) {

    auto e = beg;
    
    // Case 1: random access iterator
    if constexpr(std::is_same_v<category, std::random_access_iterator_tag>) {
      size_t r = std::distance(beg, end);
      std::advance(e, std::min(r, g));
    }
    // Case 2: non-random access iterator
    else {
      for(size_t i=0; i<g && e != end; ++e, ++i);
    }
      
    // Create a task
    auto [task, future] = emplace([beg, e, bop, uop] () mutable {
      auto init = uop(*beg);
      for(++beg; beg != e; ++beg) {
        init = bop(std::move(init), uop(*beg));          
      }
      return init;
    });
    source.precede(task);
    task.precede(target);
    futures.push_back(std::move(future));

    // adjust the pointer
    beg = e;
  }

  // target synchronizer
  target.work([&result, futures=MoC{std::move(futures)}, bop] () {
    for(auto& fu : futures.object) {
      result = bop(std::move(result), fu.get());
    }
  });

  return std::make_pair(source, target); 
}

// Function: transform_reduce    
template <typename Node>
template <typename I, typename T, typename B, typename P, typename U>
auto BasicFlowBuilder<Node>::transform_reduce(
  I beg, I end, T& result, B&& bop, P&& pop, U&& uop
) {

  using category = typename std::iterator_traits<I>::iterator_category;
  
  // Even partition
  size_t d = std::distance(beg, end);
  size_t w = std::max(size_t{1}, _num_workers);
  size_t g = std::max((d + w - 1) / w, size_t{2});

  auto source = placeholder();
  auto target = placeholder();

  std::vector<std::future<T>> futures;

  while(beg != end) {

    auto e = beg;
    
    // Case 1: random access iterator
    if constexpr(std::is_same_v<category, std::random_access_iterator_tag>) {
      size_t r = std::distance(beg, end);
      std::advance(e, std::min(r, g));
    }
    // Case 2: non-random access iterator
    else {
      for(size_t i=0; i<g && e != end; ++e, ++i);
    }
      
    // Create a task
    auto [task, future] = emplace([beg, e, uop, pop] () mutable {
      auto init = uop(*beg);
      for(++beg; beg != e; ++beg) {
        init = pop(std::move(init), *beg);
      }
      return init;
    });
    source.precede(task);
    task.precede(target);
    futures.push_back(std::move(future));

    // adjust the pointer
    beg = e;
  }

  // target synchronizer
  target.work([&result, futures=MoC{std::move(futures)}, bop] () {
    for(auto& fu : futures.object) {
      result = bop(std::move(result), fu.get());
    }
  });

  return std::make_pair(source, target); 
}


// Procedure: _linearize
template <typename Node>
template <typename L>
void BasicFlowBuilder<Node>::_linearize(L& keys) {
  (void) std::adjacent_find(
    keys.begin(), keys.end(), 
    [] (auto& from, auto& to) {
      from._node->precede(*(to._node));
      return false;
    }
  );
}

// Procedure: linearize
template <typename Node>
void BasicFlowBuilder<Node>::linearize(std::vector<Task>& keys) {
  _linearize(keys); 
}

// Procedure: linearize
template <typename Node>
void BasicFlowBuilder<Node>::linearize(std::initializer_list<Task> keys) {
  _linearize(keys);
}

// Proceduer: reduce
template <typename Node>
template <typename I, typename T, typename B>
auto BasicFlowBuilder<Node>::reduce(I beg, I end, T& result, B&& op) {
  
  using category = typename std::iterator_traits<I>::iterator_category;
  
  size_t d = std::distance(beg, end);
  size_t w = std::max(size_t{1}, _num_workers);
  size_t g = std::max((d + w - 1) / w, size_t{2});

  auto source = placeholder();
  auto target = placeholder();

  std::vector<std::future<T>> futures;
  
  while(beg != end) {

    auto e = beg;
    
    // Case 1: random access iterator
    if constexpr(std::is_same_v<category, std::random_access_iterator_tag>) {
      size_t r = std::distance(beg, end);
      std::advance(e, std::min(r, g));
    }
    // Case 2: non-random access iterator
    else {
      for(size_t i=0; i<g && e != end; ++e, ++i);
    }
      
    // Create a task
    auto [task, future] = emplace([beg, e, op] () mutable {
      auto init = *beg;
      for(++beg; beg != e; ++beg) {
        init = op(std::move(init), *beg);          
      }
      return init;
    });
    source.precede(task);
    task.precede(target);
    futures.push_back(std::move(future));

    // adjust the pointer
    beg = e;
  }
  
  // target synchronizer
  target.work([&result, futures=MoC{std::move(futures)}, op] () {
    for(auto& fu : futures.object) {
      result = op(std::move(result), fu.get());
    }
  });

  return std::make_pair(source, target); 
}

// ----------------------------------------------------------------------------

// Class: BasicSubflowBuilder
template <typename Node>
class BasicSubflowBuilder : public BasicFlowBuilder<Node> {

  using Base = BasicFlowBuilder<Node>;

  public:
    
    template <typename... Args>
    BasicSubflowBuilder(Args&&...);

    void join();
    void detach();

    bool detached() const;
    bool joined() const;

  private:

    bool _detached {false};
};

// Constructor
template <typename Node>
template <typename... Args>
BasicSubflowBuilder<Node>::BasicSubflowBuilder(Args&&... args) :
  Base {std::forward<Args>(args)...} {
}

// Procedure: join
template <typename Node>
void BasicSubflowBuilder<Node>::join() {
  _detached = false;
}

// Procedure: detach
template <typename Node>
void BasicSubflowBuilder<Node>::detach() {
  _detached = true;
}

// Function: detached
template <typename Node>
bool BasicSubflowBuilder<Node>::detached() const {
  return _detached;
}

// Function: joined
template <typename Node>
bool BasicSubflowBuilder<Node>::joined() const {
  return !_detached;
}

// ============================================================================
// Definition starts BasicTaskflow
// ============================================================================

// Class: BasicTaskflow
template <template <typename...> typename F, typename P>
class BasicTaskflow {
  
  public:

  using Node           = BasicNode<F>;
  using Work           = typename Node::Work;
  using Subwork        = typename Node::Subwork;

  using Task           = BasicTask<Node>;
  using FlowBuilder    = BasicFlowBuilder<Node>;
  using SubflowBuilder = BasicSubflowBuilder<Node>;
  using Topology       = BasicTopology<Node>;
 
    BasicTaskflow();
    BasicTaskflow(unsigned);
    ~BasicTaskflow();
    
    template <typename C>
    auto emplace(C&&);

    template <typename... C, std::enable_if_t<(sizeof...(C)>1), void>* = nullptr>
    auto emplace(C&&...);

    template <typename C>
    auto silent_emplace(C&&);

    template <typename... C, std::enable_if_t<(sizeof...(C)>1), void>* = nullptr>
    auto silent_emplace(C&&...);

    template <typename I, typename C>
    auto parallel_for(I, I, C&&, size_t = 0);

    template <typename T, typename C, std::enable_if_t<is_iterable_v<T>, void>* = nullptr>
    auto parallel_for(T&, C&&, size_t = 0);

    template <typename I, typename T, typename B>
    auto reduce(I, I, T&, B&&);

    template <typename I, typename T>
    auto reduce_min(I, I, T&);
    
    template <typename I, typename T>
    auto reduce_max(I, I, T&);

    template <typename I, typename T, typename B, typename U>
    auto transform_reduce(I, I, T&, B&&, U&&);

    template <typename I, typename T, typename B, typename R, typename U>
    auto transform_reduce(I, I, T&, B&&, R&&, U&&);

    auto placeholder();

    std::shared_future<void> dispatch();

    void precede(Task, Task);
    void linearize(std::vector<Task>&);
    void linearize(std::initializer_list<Task>);
    void broadcast(Task, std::vector<Task>&);
    void broadcast(Task, std::initializer_list<Task>);
    void gather(std::vector<Task>&, Task);
    void gather(std::initializer_list<Task>, Task);  
    void silent_dispatch();
    void wait_for_all();
    void wait_for_topologies();
    void num_workers(size_t);

    size_t num_nodes() const;
    size_t num_workers() const;
    size_t num_topologies() const;

    std::string dump() const;
    std::string dump_topologies() const;

  private:

    P _threadpool;

    std::forward_list<Node> _nodes;
    std::forward_list<Topology> _topologies;

    void _schedule(Node&);
};

// Constructor
template <template <typename...> typename F, typename P>
BasicTaskflow<F, P>::BasicTaskflow() : 
  _threadpool{std::thread::hardware_concurrency()} {
}

// Constructor
template <template <typename...> typename F, typename P>
BasicTaskflow<F, P>::BasicTaskflow(unsigned N) : _threadpool{N} {
}

// Destructor
template <template <typename...> typename F, typename P>
BasicTaskflow<F, P>::~BasicTaskflow() {
  wait_for_topologies();
}

// Procedure: num_workers
template <template <typename...> typename F, typename P>
void BasicTaskflow<F, P>::num_workers(size_t W) {
  _threadpool.shutdown();
  _threadpool.spawn(W);
}

// Function: num_nodes
template <template <typename...> typename F, typename P>
size_t BasicTaskflow<F, P>::num_nodes() const {
  //return _nodes.size();
  return std::distance(_nodes.begin(), _nodes.end());
}

// Function: num_workers
template <template <typename...> typename F, typename P>
size_t BasicTaskflow<F, P>::num_workers() const {
  return _threadpool.num_workers();
}

// Function: num_topologies
template <template <typename...> typename F, typename P>
size_t BasicTaskflow<F, P>::num_topologies() const {
  return _topologies.size();
}

// Procedure: precede
template <template <typename...> typename F, typename P>
void BasicTaskflow<F, P>::precede(Task from, Task to) {
  from._node->precede(*(to._node));
}

// Procedure: linearize
template <template <typename...> typename F, typename P>
void BasicTaskflow<F, P>::linearize(std::vector<Task>& keys) {
  FlowBuilder(_nodes, num_workers()).linearize(keys);
}

// Procedure: linearize
template <template <typename...> typename F, typename P>
void BasicTaskflow<F, P>::linearize(std::initializer_list<Task> keys) {
  FlowBuilder(_nodes, num_workers()).linearize(keys);
}

// Procedure: broadcast
template <template <typename...> typename F, typename P>
void BasicTaskflow<F, P>::broadcast(Task from, std::vector<Task>& keys) {
  from.broadcast(keys);
}

// Procedure: broadcast
template <template <typename...> typename F, typename P>
void BasicTaskflow<F, P>::broadcast(Task from, std::initializer_list<Task> keys) {
  from.broadcast(keys);
}


// Function: gather
template <template <typename...> typename F, typename P>
void BasicTaskflow<F, P>::gather(std::vector<Task>& keys, Task to) {
  to.gather(keys);
}

// Function: gather
template <template <typename...> typename F, typename P>
void BasicTaskflow<F, P>::gather(std::initializer_list<Task> keys, Task to) {
  to.gather(keys);
}

// Procedure: silent_dispatch 
template <template <typename...> typename F, typename P>
void BasicTaskflow<F, P>::silent_dispatch() {

  if(_nodes.empty()) return;

  auto& topology = _topologies.emplace_front(std::move(_nodes));

  // Start the taskflow
  _schedule(topology._source);
}

// Procedure: dispatch 
template <template <typename...> typename F, typename P>
std::shared_future<void> BasicTaskflow<F, P>::dispatch() {

  if(_nodes.empty()) {
    return std::async(std::launch::deferred, [](){}).share();
  }

  auto& topology = _topologies.emplace_front(std::move(_nodes));

  // Start the taskflow
  _schedule(topology._source);
  
  return topology._future;
}

// Procedure: wait_for_all
template <template <typename...> typename F, typename P>
void BasicTaskflow<F, P>::wait_for_all() {
  if(!_nodes.empty()) {
    silent_dispatch();
  }
  wait_for_topologies();
}

// Procedure: wait_for_topologies
template <template <typename...> typename F, typename P>
void BasicTaskflow<F, P>::wait_for_topologies() {
  for(auto& t: _topologies){
    t._future.get();
  }
  _topologies.clear();
}

// Function: placeholder
template <template <typename...> typename F, typename P>
auto BasicTaskflow<F, P>::placeholder() {
  return FlowBuilder(_nodes, num_workers()).placeholder();
}

// Function: silent_emplace
template <template <typename...> typename F, typename P>
template <typename C>
auto BasicTaskflow<F, P>::silent_emplace(C&& c) {
  return FlowBuilder(_nodes, num_workers()).silent_emplace(std::forward<C>(c));
}

// Function: silent_emplace
template <template <typename...> typename F, typename P>
template <typename... C, std::enable_if_t<(sizeof...(C)>1), void>*>
auto BasicTaskflow<F, P>::silent_emplace(C&&... cs) {
  return FlowBuilder(_nodes, num_workers()).silent_emplace(std::forward<C>(cs)...);
}

// Function: emplace
template <template <typename...> typename F, typename P>
template <typename C>
auto BasicTaskflow<F, P>::emplace(C&& c) {
  return FlowBuilder(_nodes, num_workers()).emplace(std::forward<C>(c));
}

// Function: emplace
template <template <typename...> typename F, typename P>
template <typename... C, std::enable_if_t<(sizeof...(C)>1), void>*>
auto BasicTaskflow<F, P>::emplace(C&&... cs) {
  return FlowBuilder(_nodes, num_workers()).emplace(std::forward<C>(cs)...);
}

// Function: parallel_for    
template <template <typename...> typename F, typename P>
template <typename I, typename C>
auto BasicTaskflow<F, P>::parallel_for(I beg, I end, C&& c, size_t g) {
  return FlowBuilder(_nodes, num_workers()).parallel_for(
    beg, end, std::forward<C>(c), g
  );
}

// Function: parallel_for
template <template <typename...> typename F, typename P>
template <typename T, typename C, std::enable_if_t<is_iterable_v<T>, void>*>
auto BasicTaskflow<F, P>::parallel_for(T& t, C&& c, size_t group) {
  return FlowBuilder(_nodes, num_workers()).parallel_for(
    t, std::forward<C>(c), group
  );
}

// Function: reduce 
template <template <typename...> typename F, typename P>
template <typename I, typename T, typename B>
auto BasicTaskflow<F, P>::reduce(I beg, I end, T& result, B&& op) {
  return FlowBuilder(_nodes, num_workers()).reduce(
    beg, end, result, std::forward<B>(op)
  );
}

// Function: reduce_min
// Find the minimum element over a range of items.
template <template <typename...> typename F, typename P>
template <typename I, typename T>
auto BasicTaskflow<F, P>::reduce_min(I beg, I end, T& result) {
  return reduce(beg, end, result, [] (const auto& l, const auto& r) {
    return std::min(l, r);
  });
}

// Function: reduce_max
// Find the maximum element over a range of items.
template <template <typename...> typename F, typename P>
template <typename I, typename T>
auto BasicTaskflow<F, P>::reduce_max(I beg, I end, T& result) {
  return reduce(beg, end, result, [] (const auto& l, const auto& r) {
    return std::max(l, r);
  });
}

// Function: transform_reduce    
template <template <typename...> typename F, typename P>
template <typename I, typename T, typename B, typename U>
auto BasicTaskflow<F, P>::transform_reduce(
  I beg, I end, T& result, B&& bop, U&& uop
) {
  return FlowBuilder(_nodes, num_workers()).transform_reduce(
    beg, end, result, std::forward<B>(bop), std::forward<U>(uop)
  );
}

// Function: transform_reduce    
template <template <typename...> typename F, typename P>
template <typename I, typename T, typename B, typename R, typename U>
auto BasicTaskflow<F, P>::transform_reduce(
  I beg, I end, T& result, B&& bop, R&& pop, U&& uop
) {
  return FlowBuilder(_nodes, num_workers()).transform_reduce(
    beg, end, result, std::forward<B>(bop), std::forward<R>(pop), std::forward<U>(uop)
  );
}

// Procedure: _schedule
// The main procedure to schedule a give task node.
// Each task node has two types of tasks - regular and subflow.
template <template <typename...> typename F, typename P>
void BasicTaskflow<F, P>::_schedule(Node& node) {

  _threadpool.silent_async([this, &node](){

    // Here we need to fetch the num_successors first to avoid the invalid memory
    // access caused by topology clear.
    const auto num_successors = node.num_successors();
    
    // regular node type
    // The default node work type. We only need to execute the callback if any.
    if(auto index=node._work.index(); index == 0) {
      if(auto &f = std::get<Work>(node._work); f != nullptr){
        std::invoke(f);
      }
    }
    // subflow node type 
    // The first time we enter into the subflow context, "subnodes" must be empty.
    // After executing the user's callback on subflow, there will be at least one
    // node node used as "super source". The second time we enter this context we 
    // don't have to reexecute the work again.
    else {
      assert(std::holds_alternative<Subwork>(node._work));
      
      SubflowBuilder fb(node._children, num_workers());

      bool empty_graph = node._children.empty();

      std::invoke(std::get<Subwork>(node._work), fb);
      
      // Need to create a subflow
      if(empty_graph) {

        auto& S = node._children.emplace_front([](){});

        S._topology = node._topology;

        for(auto i = std::next(node._children.begin()); i != node._children.end(); ++i) {

          i->_topology = node._topology;

          if(i->num_successors() == 0) {
            i->precede(fb.detached() ? node._topology->_target : node);
          }

          if(i->num_dependents() == 0) {
            S.precede(*i);
          }
        }
        
        // this is for the case where subflow graph might be empty
        if(!fb.detached()) {
          S.precede(node);
        }

        _schedule(S);

        if(!fb.detached()) {
          return;
        }
      }
    }
    
    // At this point, the node/node storage might be destructed.
    for(size_t i=0; i<num_successors; ++i) {
      if(--(node._successors[i]->_dependents) == 0) {
        _schedule(*(node._successors[i]));
      }
    }
  });
}

// Function: dump_topology
template <template <typename...> typename F, typename P>
std::string BasicTaskflow<F, P>::dump_topologies() const {
  
  std::ostringstream os;

  for(const auto& tpg : _topologies) {
    tpg._dump(os);
  }
  
  return os.str();
}

// Function: dump
// Dumps the taskflow in graphviz. The result can be viewed at http://www.webgraphviz.com/.
template <template <typename...> typename F, typename P>
std::string BasicTaskflow<F, P>::dump() const {

  std::ostringstream os;

  os << "digraph Taskflow {\n";
  
  for(const auto& node : _nodes) {
    node._dump(os);
  }

  os << "}\n";
  
  return os.str();
}

//-----------------------------------------------------------------------------

using Taskflow = BasicTaskflow<std::function, SpeculativeThreadpool>;
using Task = typename Taskflow::Task;
using SubflowBuilder = typename Taskflow::SubflowBuilder;

};  // end of namespace tf. ---------------------------------------------------





