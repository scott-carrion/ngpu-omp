// created 10/04/2019 by Brennan Young
// modified 10/08/2019 by Brennan Young
// - default constructor now takes an optional argument to set value.

#ifndef YOUNG_GRIDCELL_20191004
#define YOUNG_GRIDCELL_20191004

template <typename T>
class GridCell {
private:
public:
    T value;
    int i, j;
    
    // default constructor
    GridCell ( int row=-1, int col=-1, T val=0.0 )
        : i(row), j(col), value(val)
        {}
    
    // copy constructor
    template <typename U>
    GridCell ( const GridCell<U> & cell )
        : i(cell.i), j(cell.j), value((T) cell.value)
        {}
    
    // destructor
    ~GridCell () {}
    
    // assignment operator
    template <typename U>
    GridCell<T> & operator= ( const GridCell<U> & cell )
    {
        i = cell.i;
        j = cell.j;
        value = (T) cell.value;
    }
    
    // logical operators
    template <typename U>
    bool operator< ( const GridCell<U> & cell ) const
    {
        if ( value < cell.value ) return true;
        if ( cell.value < value ) return false;
        if ( i < cell.i ) return true;
        if ( cell.i < i ) return false;
        return j < cell.j;
    }
};

#endif // YOUNG_GRIDCELL_20191004