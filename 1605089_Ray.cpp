
#include "1605089_point.cpp"

class Ray {

    public:
        struct point start;
        struct point dir;

        Ray(struct point s, struct point d);
        ~Ray(){};

};

Ray::Ray(struct point s, struct point d){

    start.x = s.x;
    start.y = s.y;
    start.z = s.z;

    dir.x = d.x;
    dir.y = d.y;
    dir.z = d.z;

}