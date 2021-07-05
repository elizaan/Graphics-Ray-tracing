#include "1605089_Ray.cpp"

class Light{

    public:

        struct point position;
        double colors[3];

        Light(struct point pos, double c1, double c2, double c3);

		~Light(){};
        
        

};

Light ::Light(struct point pos, double c1, double c2, double c3){

    position.x = pos.x;
    position.y = pos.y;
    position.z = pos.z;
    colors[0] = c1;
    colors[1] = c2;
    colors[2] = c3;

}