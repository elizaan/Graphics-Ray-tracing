#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <initializer_list>
#include "bitmap_image.hpp"
#include "1605089_Light.cpp"


// #include <GL/glext.h>

#define pi (2*acos(0.0))
using namespace std;



void drawSphere(double radius,int slices,int stacks, double colr, double colg, double colb)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	// int color = 0;
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
        
		
		for(j=0;j<slices;j++)
		{
			glColor3f(colr, colg, colb);
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();

			// color = 1- color;
		}
		
	}
}

void drawTriangle(struct point a, struct point b, struct point c,double colr, double colg, double colb){
    glBegin(GL_TRIANGLES);
        {
            glColor3d(colr, colg, colb);
            glVertex3d(a.x, a.y, a.z);
            glVertex3d(b.x, b.y, b.z);
            glVertex3d(c.x, c.y, c.z);
        }
        glEnd();
}

class Object {

    public:


        struct point reference_point; // should have x, y, z
        double height, width, length;
        double col[3];
        double coEfficients[4];
    // reflection coefficients
        int shine; // exponent term of specular component
        Object();
   
        virtual void draw() = 0;
        virtual void printA() = 0;

        virtual double intersect(Ray r, double *Color, int level){
            return -1.0;
        }
        void setColor(double r, double g, double b);
        void setShine(int s);
        void setCoEfficients(double ambient, double gdiffuse, double specular, double recursive_reflection_coefficient);

        virtual ~Object() {}
};

Object:: Object(){


}

void Object:: setColor(double r, double g, double b){
    col[0] = r;
    col[1] = g;
    col[2] = b;
}

void Object:: setShine(int s){
    shine = s;
}

void Object:: setCoEfficients(double ambient, double gdiffuse, double specular, double recursive_reflection_coefficient){

    coEfficients[0] = ambient;
    coEfficients[1] = gdiffuse;
    coEfficients[2] = specular;
    coEfficients[3] = recursive_reflection_coefficient;

}


extern vector <Light> lights;
extern vector <Object*> objects;
extern int levelOfRecursion;

class Sphere: public Object{
    public:

        // point center;

        Sphere(point center, double radius);
        void draw();
        void printA();
        double intersect(Ray r, double *Color, int level);
        
};

Sphere::Sphere(struct point center, double radius){


    reference_point.x = center.x;
    reference_point.y = center.y;
    reference_point.z = center.z;
    length = radius;


}

void Sphere::draw(){

    glPushMatrix();

    glTranslatef(reference_point.x, reference_point.y, reference_point.z);

    drawSphere(length,50,50, col[0],col[1], col[2]);
    
    glPopMatrix();

    

}

void Sphere:: printA(){
    cout << " "<<endl;
}

double Sphere:: intersect(Ray r, double *Color, int level){

    struct point R0, Rd, intersectionPoint;
    R0.x = r.start.x - reference_point.x;
    R0.y = r.start.y - reference_point.y;
    R0.z = r.start.z - reference_point.z;

    Rd.x = r.dir.x;
    Rd.y = r.dir.y;
    Rd.z = r.dir.z;

    double a, b, c, D, t1, t2, T;

    a = 1;
    b = 2 * dotProduct(Rd, R0);
    c = dotProduct(R0, R0) - (length * length);

    D = (b * b) - (4 * a * c);

    if (D < 0)
    {
        return -1;
    }
    t1 = (-b + sqrt(D)) / 2 * a;
    t2 = (-b - sqrt(D)) / 2 * a;

    if (t1 < 0 && t2 < 0)
    {
        return -1;
    }
    else if (t1 > 0 && t2 > 0)
    {
        T = min(t1, t2);
    }
    else if (t1 > 0 && t2 < 0)
    {
        T = t1;
    }
    else if (t1 < 0 && t2 > 0)
    {
        T = t2;
    }

    if (level == 0)
    {
        return T;
    }

    intersectionPoint.x = r.start.x + T * r.dir.x;
    intersectionPoint.y = r.start.y + T * r.dir.y;
    intersectionPoint.z = r.start.z + T * r.dir.z;

    for (int i = 0; i < 3; i++)
    {
        Color[i] = col[i] * coEfficients[0];
    }

    //illumination starts here

    struct point normal;
    normal.x = intersectionPoint.x - reference_point.x; // from reference point to intersection point
    normal.y = intersectionPoint.y - reference_point.y;
    normal.z = intersectionPoint.z - reference_point.z;

    normal = normalize(normal);
    double distance = sqrt((r.start.x - reference_point.x) *(r.start.x - reference_point.x) +
                    (r.start.y - reference_point.y) *(r.start.y - reference_point.y)+ 
                    (r.start.z - reference_point.z)* (r.start.z - reference_point.z));

    if(distance < length){ // checking if the  eye is inside of the sphere or not
        normal = negateVector(normal); // chamging the dorection
    }

    struct point reflection;

    for (int i = 0; i < lights.size(); i++)
    {
        double llength;

        struct point ldirection = {lights[i].position.x - intersectionPoint.x, lights[i].position.y - intersectionPoint.y, lights[i].position.z - intersectionPoint.z}; 
        // from intersection point to light position

        llength = sqrt(dotProduct(ldirection, ldirection));
        ldirection = normalize(ldirection);
        struct point rstart;
        rstart = {intersectionPoint.x + ldirection.x * 1, intersectionPoint.y + ldirection.y * 1, intersectionPoint.z + ldirection.z * 1};

        Ray L(rstart, ldirection);

        bool isObscured = false;

        for (int j = 0; j < objects.size(); j++)
        {

            double t = objects[j]->intersect(L, Color, 0);

            if (t > 0 && t < llength)
            {
                isObscured = true;
                break;
            }
        }

        if (isObscured == false)
        {
            double lambart = dotProduct(normal, L.dir);

            reflection = {L.dir.x - 2 * dotProduct(L.dir, normal) * normal.x, L.dir.y - 2 * dotProduct(L.dir, normal) * normal.y, L.dir.z - 2 * dotProduct(L.dir, normal) * normal.z};
            reflection = normalize(reflection);
            double rr = dotProduct(reflection, r.dir);
            double phong = pow(max(0.0, rr), shine);

            lambart = max(lambart, 0.0);

            for (int k = 0; k < 3; k++)
            {
                Color[k] = Color[k] + col[k] * lights[i].colors[k] * coEfficients[1] * lambart; // MULTIPLYING WITH OWN COLOR
                Color[k] = Color[k] + col[k] * lights[i].colors[k] * coEfficients[2] * phong;
            }

            // for (int k = 0; k < 3; k++) //color clipping
            // {
            //     if (Color[k] < 0.0)
            //     {
            //         Color[k] = 0.0;
            //     }
            //     else if (Color[k] > 1.0)
            //     {
            //         Color[k] = 1.0;
            //     }
            // }
        }
    }

    
    // Reflection starts here 

    if (level>= levelOfRecursion){
        return T;
    }

    struct point reflection2 = {r.dir.x - 2* dotProduct(r.dir, normal)*normal.x, r.dir.y - 2* dotProduct(r.dir, normal)*normal.y, r.dir.z - 2* dotProduct(r.dir, normal)*normal.z};
    reflection2 = normalize(reflection2);

    struct point startr = {intersectionPoint.x + reflection2.x * 1, intersectionPoint.y + reflection2.y *1, intersectionPoint.z+reflection2.z*1};

    Ray reflect(startr, reflection2);

    double tc;
	double tMin = 999999;
	int nearest = -1;

    double *reflectColor = new double[3];


    for(int k = 0; k < objects.size(); k++){

		tc = objects[k]->intersect(reflect, reflectColor, 0);
				// if( k == 7){
				// 	cout << " for "<< k <<"th object t is " << t << endl;
				// }
				
		if(tc>0 && tc< tMin){
					nearest = k;
					tMin = tc;
					
		}


	}
	if( nearest!= -1){

		double temp;
		temp = objects[nearest]->intersect(reflect, reflectColor, level+1);

        for( int l = 0; l <3 ; l++){
            Color[l]= Color[l]+  reflectColor[l] * coEfficients[3];
        }
		
        for (int l = 0; l < 3; l++){ //color clipping
            if(Color[l]<0.0){
                Color[l] = 0.0;
            }
            else if (Color[l]>1.0){
                Color[l] = 1.0;
            }
        }
				
				

	}

    return T;
}

class Triangle: public Object{
    public:
        struct point A, B, C;

        Triangle(struct point A, struct point B, struct point C){
            this->A = A;
            this->B = B;
            this->C = C;
        }

        void draw();
        void printA();
        double intersect(Ray r, double *Color, int level);

};

void Triangle::draw(){

    drawTriangle(this->A,this->B,this->C,col[0],col[1],col[2]);

}

void Triangle:: printA(){
    cout << this->A.x << " " << this->B.x << " " << this->C.x << " "<< endl;
}

double Triangle:: intersect(Ray r, double *Color, int level){

    double T;
    // from the given link
    struct point intersectionPoint;


    const double EPSILON = 0.0000001;
    struct point vertex0 = this->A;
    struct point vertex1 = this->B;  
    struct point vertex2 = this->C;
    struct point edge1, edge2, h, s, q;
    double a,f,u,v;
    edge1.x = vertex1.x - vertex0.x;
    edge1.y = vertex1.y - vertex0.y;
    edge1.z = vertex1.z - vertex0.z;

    edge2.x = vertex2.x - vertex0.x;
    edge2.y = vertex2.y - vertex0.y;
    edge2.z = vertex2.z - vertex0.z;

    h = crossProduct(r.dir, edge2);
    a = dotProduct(edge1,h);

    if (a > -EPSILON && a < EPSILON)
        return -1;    // This ray is parallel to this triangle.
    f = 1.0/a;
    s.x = r.start.x- vertex0.x;
    s.y = r.start.y - vertex0.y;
    s.z = r.start.z- vertex0.z;

    u = f * dotProduct(s,h);

    if (u < 0.0 || u > 1.0)
        return -1;
    q = crossProduct(s,edge1);
    v = f * dotProduct(r.dir,q);

    if (v < 0.0 || (u + v) > 1.0)
        return -1;
    // At this stage we can compute t to find out where the intersection point is on the line.
    T = f * dotProduct(edge2,q);
    if(level == 0){

        if (T > EPSILON) // ray intersection
        {
            intersectionPoint.x = r.start.x + r.dir.x * T;
            intersectionPoint.y = r.start.y + r.dir.y * T;
            intersectionPoint.z = r.start.z + r.dir.z * T;

            return T;
        }
        else {
            return -1;
        }            // This means that there is a line intersection but not a ray intersection.
    }
        


    for( int i = 0; i < 3; i++){
            Color[i] = col[i]* coEfficients[0];
    }

    // Illumination starts here 

    struct point normal;
    normal = crossProduct(edge1, edge2);
    
    normal = normalize(normal);

    if(dotProduct(r.dir, normal)> 0.0){
        normal = negateVector(normal);
    }

    struct point reflection;
       

    for(int i = 0; i < lights.size(); i++){
        double llength;

        struct point ldirection = {lights[i].position.x - intersectionPoint.x, lights[i].position.y - intersectionPoint.y, lights[i].position.z - intersectionPoint.z};
           
        llength = sqrt(dotProduct(ldirection, ldirection));
        ldirection = normalize(ldirection);
        struct point rstart;
        rstart = {intersectionPoint.x + ldirection.x *1, intersectionPoint.y + ldirection.y *1, intersectionPoint.z + ldirection.z *1};

        Ray L(rstart, ldirection);

        bool isObscured = false;


        for ( int j = 0; j < objects.size(); j++){

            double t = objects[j]->intersect(L, Color, 0);

            if (t>0 && t< llength){
                    isObscured = true;
                    break;
            }
        }

        if (isObscured == false){
            double lambart = dotProduct(normal, L.dir);

            reflection = {L.dir.x - 2* dotProduct(L.dir, normal)*normal.x, L.dir.y - 2* dotProduct(L.dir, normal)*normal.y, L.dir.z - 2* dotProduct(L.dir, normal)*normal.z};
            reflection = normalize(reflection);
            double rr = dotProduct(reflection, r.dir);
            double phong = pow(max(0.0,rr), shine);

            lambart = max(lambart, 0.0);

            for( int w = 0; w < 3; w++){
                Color[w] = Color[w]+ col[w]* lights[i].colors[w]* coEfficients[1]*lambart; // MULTIPLYING WITH OWN COLOR
                Color[w] = Color[w]+ col[w]* lights[w].colors[w]* coEfficients[2]*phong;

                    

            }

            // for (int w = 0; w < 3; w++){
            //     if(Color[w]<0.0){
            //         Color[w] = 0.0;
            //     }
            //     else if (Color[w]>1.0){
            //         Color[w] = 1.0;
            //     }
            // }



        }


    }


    // Reflection starts here 

    if (level>= levelOfRecursion){
        return T;
    }

    struct point reflection2 = {r.dir.x - 2* dotProduct(r.dir, normal)*normal.x, r.dir.y - 2* dotProduct(r.dir, normal)*normal.y, r.dir.z - 2* dotProduct(r.dir, normal)*normal.z};
    reflection2 = normalize(reflection2);

    struct point startr = {intersectionPoint.x + reflection2.x * 1, intersectionPoint.y + reflection2.y *1, intersectionPoint.z+reflection2.z*1};

    Ray reflect(startr, reflection2);

    double tc;
	double tMin = 999999;
	int nearest = -1;

    double *reflectColor = new double[3];


    for(int k = 0; k < objects.size(); k++){

		tc = objects[k]->intersect(reflect, reflectColor, 0);
				// if( k == 7){
				// 	cout << " for "<< k <<"th object t is " << t << endl;
				// }
				
		if(tc>0 && tc< tMin){
					nearest = k;
					tMin = tc;
					
		}


	}
	if( nearest!= -1){

		double temp;
		temp = objects[nearest]->intersect(reflect, reflectColor, level+1);

        for( int l = 0; l <3 ; l++){
            Color[l]= Color[l]+  reflectColor[l] * coEfficients[3];
        }
		
        for (int l = 0; l < 3; l++){
            if(Color[l]<0.0){
                Color[l] = 0.0;
            }
            else if (Color[l]>1.0){
                Color[l] = 1.0;
            }
        }
				
				

	}

    return T;
}

class Floor: public Object{

    public:

    int n;
    // struct point reference;
    double flength, tlength;
    
    

    Floor(double fl, double tl);
    void draw();
    double intersect(Ray r, double *Color, int level);
    void printA();
};

Floor::Floor(double fl, double tl){

        flength = fl;
        tlength = tl;

        reference_point.x = -flength/2.0;
        reference_point.y = -flength/2.0;
        reference_point.z = 0.0;
        n = flength/tlength;

        
}

void Floor::draw(){

    glBegin(GL_QUADS);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((i + j) % 2 == 0){
                    glColor3d(1.0, 1.0, 1.0);
                }
                    
                else{

                    glColor3d(0.0, 0.0, 0.0);

                }
                    
                glVertex3d(reference_point.x + tlength * i, reference_point.y + tlength * j, reference_point.z);
                glVertex3d(reference_point.x + tlength* (i + 1), reference_point.y + tlength* j, reference_point.z);
                glVertex3d(reference_point.x + tlength * (i + 1), reference_point.y + tlength * (j + 1), reference_point.z);
                glVertex3d(reference_point.x + tlength * i, reference_point.y + tlength * (j + 1), reference_point.z);
            }
        }
        glEnd();

}

double Floor:: intersect(Ray r, double *Color, int level){


    double T;
    double D = 0;
    struct point normal, intersectionPoint;
    int i1, j1;
    normal.x = 0.0;
    normal.y =0.0;
    normal.z= 1.0;

    double denom = dotProduct(normal, r.dir);
    if(denom == 0){
        return -1;
    }


    T = -1.0*(dotProduct(normal,r.start)) / (dotProduct(normal,r.dir));
    // T = -(r.start.z/ r.dir.z);
    

    if(T < 0.0){
        return -1;
    }

    if(level == 0){
        return T;
    }

    intersectionPoint.x = r.start.x + T*r.dir.x;
    intersectionPoint.y = r.start.y + T*r.dir.y;
    intersectionPoint.z = r.start.z + T*r.dir.z;


    double x1 = reference_point.x;
    double y1 = reference_point.y;
    double x2 = -reference_point.x;
    double y2 = -reference_point.y; 

    if(intersectionPoint.x<x1 || intersectionPoint.x > x2 || intersectionPoint.y < y1 || intersectionPoint.y > y2){

        return -1;

    }

    i1 = (int)((intersectionPoint.x - reference_point.x)/ tlength);
    j1 = (int)((intersectionPoint.y - reference_point.y)/ tlength);

    if ((i1+j1)%2 == 0){
        this->col[0] = 1.0;
        this->col[1] = 1.0;
        this->col[2] = 1.0;
    }

    else{
        this->col[0] = 0.0;
        this->col[1] = 0.0;
        this->col[2] = 0.0;
    }

    // this->col[0] = Color[0];
    // this->col[1] = Color[1];
    // this->col[2] = Color[2];
    

    for( int k = 0; k < 3; k++){
            Color[k] = col[k]* coEfficients[0];
    }
    
    
    normal = normalize(normal);

    if(dotProduct(r.dir, normal)> 0.0){
        normal = negateVector(normal);
    }

    struct point reflection;
       

    for(int l = 0; l < lights.size(); l++){
        double llength;

        struct point ldirection = {lights[l].position.x - intersectionPoint.x, lights[l].position.y - intersectionPoint.y, lights[l].position.z - intersectionPoint.z};
           
        llength = sqrt(dotProduct(ldirection, ldirection));
        ldirection = normalize(ldirection);
        struct point rstart;
        rstart = {intersectionPoint.x + ldirection.x *1, intersectionPoint.y + ldirection.y *1, intersectionPoint.z + ldirection.z *1};

        Ray L(rstart, ldirection);

        bool isObscured = false;


        for ( int m = 0; m < objects.size(); m++){

            double t = objects[m]->intersect(L, Color, 0);

            if (t>0 && t< llength){
                    isObscured = true;
                    break;
            }
        }

        if (isObscured == false){
            double lambart = dotProduct(normal, L.dir);

            reflection = {L.dir.x - 2* dotProduct(L.dir, normal)*normal.x, L.dir.y - 2* dotProduct(L.dir, normal)*normal.y, L.dir.z - 2* dotProduct(L.dir, normal)*normal.z};
            reflection = normalize(reflection);
            double rr = dotProduct(reflection, r.dir);
            // double rr = dotProduct(reflection,negateVector(r.dir));
            double phong = pow(max(0.0,rr), shine);

            lambart = max(lambart, 0.0);

            for( int w = 0; w < 3; w++){
                Color[w] = Color[w]+ col[w]* lights[l].colors[w]* coEfficients[1]*lambart; // MULTIPLYING WITH OWN COLOR
                Color[w] = Color[w]+ col[w]* lights[l].colors[w]* coEfficients[2]*phong;

                    

            }

            // for (int w = 0; w < 3; w++){
            //     if(Color[w]<0.0){
            //         Color[w] = 0.0;
            //     }
            //     else if (Color[w]>1.0){
            //         Color[w] = 1.0;
            //     }
            // }



        }


    }

   
    // Reflection starts here 

    if (level>= levelOfRecursion){
        return T;
    }

    struct point reflection2 = {r.dir.x - 2* dotProduct(r.dir, normal)*normal.x, r.dir.y - 2* dotProduct(r.dir, normal)*normal.y, r.dir.z - 2* dotProduct(r.dir, normal)*normal.z};
    reflection2 = normalize(reflection2);

    struct point startr = {intersectionPoint.x + reflection2.x * 1, intersectionPoint.y + reflection2.y *1, intersectionPoint.z+reflection2.z*1};

    Ray reflect(startr, reflection2);

    double tc;
	double tMin = 999999;
	int nearest = -1;

    double *reflectColor = new double[3];


    for(int k = 0; k < objects.size(); k++){

		tc = objects[k]->intersect(reflect, reflectColor, 0);
				// if( k == 7){
				// 	cout << " for "<< k <<"th object t is " << t << endl;
				// }
				
		if(tc>0 && tc< tMin){
					nearest = k;
					tMin = tc;
					
		}


	}
	if( nearest!= -1){

		double temp;
		temp = objects[nearest]->intersect(reflect, reflectColor, level+1);

        for( int l = 0; l <3 ; l++){
            Color[l]= Color[l]+  reflectColor[l] * coEfficients[3];
        }
		
        for (int l = 0; l < 3; l++){
            if(Color[l]<0.0){
                Color[l] = 0.0;
            }
            else if (Color[l]>1.0){
                Color[l] = 1.0;
            }
        }
				
				

	}


    return T;




}

void Floor:: printA(){
    cout << " "<<endl;
}

class General: public Object{

    public:
        double A, B, C, D, E, F, G, H, I, J;
        General(double a, double b, double c, double d, double e, double f, double g, double h, double i, double j,struct point ref, double l, double w, double hight);
        void draw();
        bool inside(double t, Ray r);
        void printA();
        double intersect(Ray r, double *Color, int level);

};

General::General(double a, double b, double c, double d, double e, double f, double g, double h, double i, double j,struct point ref, double l, double w, double hight){

    A = a;
    B = b;
    C = c;
    D = d;
    E = e;
    F = f;
    G = g;
    H = h;
    I = i;
    J = j;
    reference_point.x = ref.x;
    reference_point.y = ref.y;
    reference_point.z = ref.z;
    length = l;
    width = w;
    height = hight;

}

void General:: printA(){

    cout<< "The coefficients  "<< A<<" " <<B <<" "<<C <<" "<<D<<" " <<E <<" "<<F <<" "<<G<<" " <<H<<" " <<I<<" "<< J<<endl;
    cout<< "colors "<<col[0]<<" "<<col[1]<<" "<<col[2]<<endl;
    cout<< "lIGHTENNG COEFFS " << coEfficients[0]<< " "<< coEfficients[1]<<" "<< coEfficients[2]<<" "<< coEfficients[3]<<endl;
    cout<<"reference points "<< reference_point.x<<" " << reference_point.y<< " "<< reference_point.z <<endl;
    cout<< " length and width and height "<< length <<" "<< width <<" "<< height<< " "<<endl;
    cout<< "Shine "<< shine<< endl;
}

void General:: draw(){

}
bool General:: inside(double t, Ray r){

    
    bool in = true; 
    struct point p;

    p.x = r.start.x+ t*r.dir.x;
    p.y = r.start.y+ t*r.dir.y;
    p.z = r.start.z+ t*r.dir.z;

    if (length != 0)
    { 
        if (abs(p.x) >= abs(reference_point.x) && abs(p.x) <= abs(reference_point.x) + length)
        {
            in = true;
        }
        else{
            return false;
        }
           
    }

    if (width != 0)
    { 
        if (abs(p.y) >= abs(reference_point.y) && abs(p.y) <= abs(reference_point.y)+ width)
        {
            in = true;
        }
        else{
            return false;
        }
            
    }

    if (height != 0)
    {
        if (abs(p.z) < abs(reference_point.z) || abs(p.z) > abs(reference_point.z) + height)
        { 
            return false;
        }
        // else
        //     in = true;
    }

    return in;

}
double General:: intersect(Ray r, double *Color, int level){

    double a,b,c,D,t1,t2,T;
    T = -1;
    bool t1Inside, t2Inside;
    struct point intersectionPoint;

    double xd = r.dir.x;
    double yd = r.dir.y;
    double zd = r.dir.z;

    double xo = r.start.x;
    double yo = r.start.y;
    double zo = r.start.z;


    a = this->A*xd*xd + this->B*yd*yd + this->C*zd*zd + this->D*xd*yd + this->E*xd*zd + this->F*yd*zd;

    b = 2*this->A*xo*xd + 2*this->B*yo*yd + 2*this->C*zo*zd + this->D*(xo*yd + yo*xd) +this->E*(xo*zd + zo*xd) + this->F*(yo*zd + yd*zo) + this->G*xd + this->H*yd + this->I*zd;

    c = this->A*xo*xo + this->B*yo*yo + this->C*zo*zo + this->D*xo*yo + this->E*xo*zo + this->F*yo*zo + this->G*xo + this->H*yo + this->I*zo + this->J;

    D = (b*b) - (4*a*c);
    
    if(a == 0.0){

        T = -c/b;

        if(T< 0.0){
            return -1;
        }

        if(inside(T, r)){
            return T;
        }
        else{
            return -1;
        }

    }

    if(D < 0.0) {
        return -1;
    }
        
    if(D == 0.0) {
        T = -b/(2*a);

        if(T< 0.0){
            return -1;
        }

        if(inside(T, r)){
            return T;
        }
        else{
            return -1;
        }
    
    }
    if(D > 0.0){
        t1 = (-b-sqrt(D))/(2*a);
        t2 = (-b+sqrt(D))/(2*a);

        if(t1 < 0.0 && t2 < 0.0){
            return -1;
        }
        else if(t1>0.0 && t2<0.0){

            // t2Inside = false;
            t1Inside = inside(t1, r);

            if(t1Inside){

                T = t1;

            }
            else{
                return -1;
            }

            
        }
        else if(t2>0.0 && t1<0.0){

            // t1Inside = false;
            t2Inside = inside(t2, r);
            if(t2Inside){

                T = t2;

            }
            else{
                return -1;
            }

            
        }

        else{
             
            t1Inside = inside(t1, r);
            t2Inside = inside(t2, r);

            if(t1Inside == false && t2Inside == false){
                return -1;
            }

            else if (t1Inside == true && t2Inside == false){
                T = t1;
            }

            else if(t1Inside == false && t2Inside == true){
                T = t2;
            }

            else{
                T = t1;
            }
            
        }
        




        

    }

    // if(T==0)
    // cout << " T : "<< T <<endl;
    if(level == 0){


        return T;
    }


    intersectionPoint.x = r.start.x + T*r.dir.x;
    intersectionPoint.y = r.start.y + T*r.dir.y;
    intersectionPoint.z = r.start.z + T*r.dir.z;



    for(int i = 0; i < 3; i++){
        Color[i] = col[i]* coEfficients[0];
    }

    //ILlumination starts here

    struct point normal;

    normal.x = 2*this->A*intersectionPoint.x + this->D*intersectionPoint.y + this->E*intersectionPoint.z+ this->G; 
    // partial differentiation
    normal.y = 2*this->B*intersectionPoint.y + this->D*intersectionPoint.x + this->F*intersectionPoint.z+ this->H;
    normal.z = 2*this->C*intersectionPoint.z + this->E*intersectionPoint.x + this->F*intersectionPoint.y+ this->I;


    normal = normalize(normal);

    if(dotProduct(r.dir, normal)> 0.0){
        normal = negateVector(normal);
    }

    struct point reflection;
       

    for(int i = 0; i < lights.size(); i++){
        double llength;

        struct point ldirection = {lights[i].position.x - intersectionPoint.x, lights[i].position.y - intersectionPoint.y, lights[i].position.z - intersectionPoint.z};
           
        llength = sqrt(dotProduct(ldirection, ldirection));
        ldirection = normalize(ldirection);
        struct point rstart;
        rstart = {intersectionPoint.x + ldirection.x *1, intersectionPoint.y + ldirection.y *1, intersectionPoint.z + ldirection.z *1};

        Ray L(rstart, ldirection);

        bool isObscured = false;


        for ( int j = 0; j < objects.size(); j++){

            double t = objects[j]->intersect(L, Color, 0);

            if (t>0 && t< llength){
                    isObscured = true;
                    break;
            }
        }

        if (isObscured == false){
            double lambart = dotProduct(normal, L.dir);

            reflection = {L.dir.x - 2* dotProduct(L.dir, normal)*normal.x, L.dir.y - 2* dotProduct(L.dir, normal)*normal.y, L.dir.z - 2* dotProduct(L.dir, normal)*normal.z};
            reflection = normalize(reflection);
            double rr = dotProduct(reflection, r.dir);
            // double rr = dotProduct(reflection,negateVector(r.dir));
            double phong = pow(max(0.0,rr), shine);

            lambart = max(lambart, 0.0);

            for( int k = 0; k < 3; k++){
                Color[k] = Color[k]+ this->col[k]* lights[i].colors[k]* coEfficients[1]*lambart; // MULTIPLYING WITH OWN COLOR
                Color[k] = Color[k]+ this->col[k]* lights[i].colors[k]* coEfficients[2]*phong;

                    

            }

            // for (int k = 0; k < 3; k++){
            //     if(Color[k]<0.0){
            //         Color[k] = 0.0;
            //     }
            //     else if (Color[k]>1.0){
            //         Color[k] = 1.0;
            //     }
            // }



        }


    }

    // Reflection starts here 

    if (level>= levelOfRecursion){
        return T;
    }

    struct point reflection2 = {r.dir.x - 2* dotProduct(r.dir, normal)*normal.x, r.dir.y - 2* dotProduct(r.dir, normal)*normal.y, r.dir.z - 2* dotProduct(r.dir, normal)*normal.z};
    reflection2 = normalize(reflection2);

    struct point startr = {intersectionPoint.x + reflection2.x * 1, intersectionPoint.y + reflection2.y *1, intersectionPoint.z+reflection2.z*1};

    Ray reflect(startr, reflection2);

    double tc;
	double tMin = 999999;
	int nearest = -1;

    double *reflectColor = new double[3];


    for(int k = 0; k < objects.size(); k++){

		tc = objects[k]->intersect(reflect, reflectColor, 0);
				// if( k == 7){
				// 	cout << " for "<< k <<"th object t is " << t << endl;
				// }
				
		if(tc>0 && tc< tMin){
					nearest = k;
					tMin = tc;
					
		}


	}
	if( nearest!= -1){

		double temp;
		temp = objects[nearest]->intersect(reflect, reflectColor, level+1);

        for( int l = 0; l <3 ; l++){
            Color[l]= Color[l]+  reflectColor[l] * coEfficients[3];
        }
		
        for (int l = 0; l < 3; l++){
            if(Color[l]<0.0){
                Color[l] = 0.0;
            }
            else if (Color[l]>1.0){
                Color[l] = 1.0;
            }
        }
				
				

	}

    return T;
}
