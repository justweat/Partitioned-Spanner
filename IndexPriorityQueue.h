//
// Created by justin on 6/12/22.
//

#ifndef INDEXPRIORITYQUEUE_INDEXPRIORITYQUEUE_H
#define INDEXPRIORITYQUEUE_INDEXPRIORITYQUEUE_H

#include <vector>
#include <unordered_map>
#include <iterator>
#include <functional>
#include <iostream>

namespace IndexPQ{
    using namespace std;

    template<class K, class V>
    function<bool(const V&, const V&)> MinHeapComparator =
            [](const V& o1, const V& o2){
                return o1 < o2;
            };

    template<class K, class V>
    function<bool(const V&, const V&)> MaxHeapComparator =
            [](const V& o1, const V& o2){
                return o1 > o2;
            };

    enum class IndexPQType {MinHeap, MaxHeap};

    template<class K, class V>
    class IndexPriorityQueue {
    public:

        IndexPriorityQueue(const vector<K> &keys, const vector<V> &vals) {

            if(keys.size() != vals.size()){
                throw logic_error("The number of keys does not match the number of values\n");
            }

            this->_heap.template emplace_back();

            size_t n = vals.size();

            for(size_t i {1}; i <= n; ++i){
                this->_heap.template emplace_back(pair<K, V>{keys[i - 1], vals[i - 1]});
            }

            for(size_t i = 1; i <= n; ++i){
                this->_keyMap.insert(pair<K, size_t>{keys[i - 1], i});
            }
        }

        IndexPriorityQueue(const vector<K> &keys,
                           const vector<V> &vals,
                           IndexPQType type)
                           : IndexPriorityQueue<K, V>(keys, vals){

            if(type == IndexPQType::MinHeap){
                this->_comparator = MinHeapComparator<K, V>;
            }else{
                this->_comparator = MaxHeapComparator<K, V>;
            }
            heapify();
        }

        IndexPriorityQueue(const vector<K> &keys,
                           const vector<V> &vals,
                           function<bool(V, V)> comparator)
                           : IndexPriorityQueue(keys, vals){
            this->_comparator = comparator;
            heapify();
        }

        void updateKey(const K& key, const V& val){

            if(!this->contains(key)){
                throw logic_error("IndexPQ does not contain key for updating");
            }
            if(this->empty()){
                throw logic_error("IndexPQ is empty");
            }

            size_t n = this->_heap.size();
            size_t pos = this->_keyMap.at(key);
            this->_heap[pos].second = val;
            size_t left_child = pos << 1;
            size_t right_child = left_child + 1;
            auto parent = pos >> 1;

            if((left_child < n && this->_comparator(this->_heap[left_child].second, this->_heap[pos].second)) || (right_child < n && this->_comparator(this->_heap[right_child].second, this->_heap[pos].second))){
                heapSink(pos);
            }else if(parent > 0 && this->_comparator(this->_heap[pos].second, this->_heap[parent].second)){
                heapSwim(pos);
            }
        }

        void push(const K &key, const V &val) {

            if(this->contains(key)){
                this->updateKey(key, val);
                return;
            }

            this->_heap.template emplace_back(pair<K,V>{key, val});
            size_t index = this->_heap.size() - 1;
            this->_keyMap.insert(pair<K, size_t>{key, index});
            heapSwim(index);
        }

        void pop(){

            if(this->empty()){
                throw logic_error("IndexPQ is empty");
            }

            popHeapMaintenance();
        }

        V frontValue() {

            if(this->empty()){
                throw logic_error("IndexPQ is empty");
            }

            return this->_heap[1].second;
        }

        K frontKey() {
            if(this->empty()){
                throw logic_error("IndexPQ is empty");
            }

            return this->_heap[1].first;
        }

        pair<K, V> frontKV() {

            if(this->empty()){
                throw logic_error("IndexPQ is empty");
            }

            return this->_heap[1];
        }

        bool contains(const K &key) {
            return this->_keyMap.find(key) != this->_keyMap.end();
        }

        bool empty() {
            return this->_heap.size() == 1;
        }

        size_t size() {
            return this->_heap.size() - 1;
        }

        V keyValue(const K &key) {

            if(!this->contains(key)){
                throw logic_error("IndexPQ does not contain key for updating");
            }

            return this->_heap[this->_keyMap.at(key)].second;
        }

        vector<K> keysWithValue(const V &val) {

            auto iter = this->_heap.begin();
            auto end = this->_heap.end();

            vector<K> keys{};
            while(iter != end){
                if(iter->second == val){
                    keys.template emplace_back(iter->second);
                }
                ++iter;
            }

            return keys;
        }

        void printIPQ() {
            cout << this->toString();
        }

        string toString() {

            IndexPriorityQueue<K, V> tempIPQ = *this;

            string ret{};
            while(!tempIPQ.empty()){
                ret += to_string(tempIPQ.frontKey()) + " -> " + to_string(tempIPQ.frontValue()) + "\n";
                tempIPQ.pop();
            }
            return ret;
        }

    private:

        void heapify() {
            size_t size = this->_heap.size();
            size_t mid = size >> 1;
            for(size_t i {mid}; i > 0; --i){
                heapSink(i);
            }
        }

        void heapSwim(size_t index) {
            size_t parent = index >> 1;
            if(parent > 0){
                if(this->_comparator(this->_heap[index].second, this->_heap[parent].second)){

                    this->_keyMap.at(this->_heap[parent].first) = index;
                    this->_keyMap.at(this->_heap[index].first) = parent;

                    swap(this->_heap[index], this->_heap[parent]);

                    heapSwim(parent);
                }
            }
        }

        void heapSink(size_t index) {
            size_t n = this->_heap.size();
            size_t left_child = index << 1;
            size_t right_child = left_child + 1;
            size_t variant;
            if(left_child < n && this->_comparator(this->_heap[left_child].second, this->_heap[index].second)){
                variant = left_child;
            }else{
                variant = index;
            }
            if(right_child < n && this->_comparator(this->_heap[right_child].second, this->_heap[variant].second)){
                variant = right_child;
            }
            if(variant != index){
                this->_keyMap.at(this->_heap[index].first) = variant;
                this->_keyMap.at(this->_heap[variant].first) = index;
                swap(this->_heap[index], this->_heap[variant]);
                heapSink(variant);
            }
        }

        void popHeapMaintenance() {

            if(this->size() == 1){
                this->_keyMap.erase(this->_heap[1].first);
                this->_heap.pop_back();
                return;
            }

            pair<K, V> back = this->_heap.back();
            this->_heap.pop_back();

            this->_keyMap.erase(this->_heap[1].first);

            this->_heap[1] = back;
            this->_keyMap.at(back.first) = 1;

            heapSink(1);
        }

        unordered_map<K, size_t> _keyMap;
        vector<pair<K, V>> _heap;
        function<bool(const V&, const V&)> _comparator;

    };
}



#endif //INDEXPRIORITYQUEUE_INDEXPRIORITYQUEUE_H
