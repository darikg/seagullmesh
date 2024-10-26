typedef boost::iterator_facade<Index_iterator<Key>, Key, std::random_access_iterator_tag, Key> Facade;

class Index_iterator : public Facade {
    public:
        Index_iterator() : i(0) {}
        Index_iterator(size_type i) : i(i) {}
    private:
        friend class boost::iterator_core_access;
        void increment() {++i;}
        void decrement() {--i;}
        void advance(size_type n) {i += n;}
        size_type distance_to(const Index_iterator& other) const {return other.n - this->n;}
        bool equal(const Index_iterator& other) const {return this->n == other.n;}
        Key dereference() const {
            // TODO need a reference to the array here, where to do the unchecked access thing?
            return hnd_;
        }

        size_type i;
};