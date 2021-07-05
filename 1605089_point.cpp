#define pi (2*acos(0.0))

struct point
{
	double x,y,z;
};

double dotProduct(struct point u, struct point v){
    return (u.x*v.x + u.y*v.y + u.z*v.z);
}

struct point crossProduct(struct point u, struct point v)
{
	struct point result;
	result.x = u.y*v.z - u.z*v.y;
	result.y = u.z*v.x - u.x*v.z;
	result.z = u.x*v.y - u.y*v.x;
 
	return result;
}

// The formula is: vector(v)*cos(angle)+ vector(vector(about)crossProduct vector(v))*sin(angle)
struct point rotation(struct point about, struct point v, double a){

	struct point scalar;
	struct point cross2;
	struct point result;

	scalar.x = v.x*(cos(a*pi/180));
	scalar.y = v.y*(cos(a*pi/180));
	scalar.z = v.z*(cos(a*pi/180));
    
	struct point cross = crossProduct(about, v); 

	cross2.x = cross.x*(sin(a*pi/180));
	cross2.y = cross.y*(sin(a*pi/180));
	cross2.z = cross.z*(sin(a*pi/180));

	result.x = scalar.x + cross2.x;
	result.y = scalar.y + cross2.y;
	result.z = scalar.z + cross2.z;
	
	return result;
}

struct point normalize(struct point u){

	double value;
	value = sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
	u.x = u.x / value;
	u.y = u.y / value;
	u.z = u.z / value;

	return u;


}

struct point negateVector(struct point u){

	struct point result;
	result.x = - u.x;
	result.y = - u.y;
	result.z = - u.z;

	return result;
}