#ifndef MODULE_H
#define MODULE_H

#include <vector>
#include <string>
#include <climits>
#include <algorithm>

class Terminal
{
public:
    Terminal(const std::string& name, int x, int y)
        : _name(name), _x1(x), _y1(y), _x2(x), _y2(y) {}
    virtual ~Terminal() {}

    const std::string& getName() const { return _name; }
    int getX1() const { return _x1; }
    int getY1() const { return _y1; }
    int getX2() const { return _x2; }
    int getY2() const { return _y2; }
    double getCenterX() const { return (_x1 + _x2) / 2.0; }
    double getCenterY() const { return (_y1 + _y2) / 2.0; }

    void setPos(int x1, int y1, int x2, int y2) {
        _x1 = x1; _y1 = y1;
        _x2 = x2; _y2 = y2;
    }

protected:
    std::string _name;
    int _x1, _y1, _x2, _y2;
};


class Block : public Terminal
{
public:
    Block(const std::string& name, int w, int h)
        : Terminal(name, 0, 0), _w(w), _h(h) {}
    ~Block() {}

    int getWidth(bool rotate = false) const { return rotate ? _h : _w; }
    int getHeight(bool rotate = false) const { return rotate ? _w : _h; }
    int getArea() const { return _w * _h; }

    void setWidth(int w) { _w = w; }
    void setHeight(int h) { _h = h; }

private:
    int _w;
    int _h;
};


class Net
{
public:
    Net() {}
    ~Net() {}

    const std::vector<Terminal*>& getTermList() const { return _termList; }
    void addTerm(Terminal* term) { _termList.push_back(term); }
    double calcHPWL() const;

private:
    std::vector<Terminal*> _termList;
};

#endif // MODULE_H
