#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>
#include <stdlib.h>

using namespace std;

struct Vec{
	double x, y, z;
	Vec() : x(0), y(0), z(0) {}
	Vec(double x, double y, double z) : x(x), y(y), z(z) {}
	Vec operator + (const Vec& v) const { return Vec(x+v.x, y+v.y, z+v.z); }
	Vec operator - (const Vec& v) const { return Vec(x-v.x, y-v.y, z-v.z); }
	Vec operator * (double d) const { return Vec(x*d, y*d, z*d); }
	Vec operator / (double d) const { return Vec(x/d, y/d, z/d); }

	Vec normalize() const{
		double dist = sqrt(x*x + y*y + z*z); 
		return Vec(x/dist, y/dist, z/dist);
	}

};

inline double dot(const Vec& v1, const Vec& v2){
	return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z); 
}

struct Ray{ // p(t) = e + td; where f(p(t)) = 0
	Vec e, d;
	Ray():e(Vec(0,0,0)), d(Vec(0,0,0)){}
	Ray(const Vec& e, const Vec& d) : e(e), d(d) {}
};

struct Sphere{
	Vec center;
	double radius;
	int reflFlag; 
	Sphere(const Vec &c, double r, int flag) : center(c), radius(r), reflFlag(flag) {}

	//unit normal
	Vec normal(const Vec& p) const{ return (p-center)/radius; }

	void setFlag(int flag) { reflFlag = flag; } 

	bool intersect(const Ray &ray, double &t) const{
		const double A = dot(ray.d, ray.d); // d . d 
		const double B = 2*dot(ray.d, ray.e - center); // 2*(d . e-c)
		const double C = dot(ray.e - center, ray.e - center) - radius*radius; // (e-c . e-c) - r*r
		double discriminant = B*B - 4*A*C; 
		//std::cout << discriminant << std::endl;

		
		if(discriminant < 0){ // IF negative
			return false; // sqrt is imaginary and no intersection
		}else{ //2 cases: 1- ray enter sphere, 2- ray leaves sphere
			discriminant = sqrt(discriminant); 
			
			const double t_minus = ((-B) - discriminant)/(2*A);
			const double t_plus = ((-B) + discriminant)/(2*A);
			
			t = (t_minus < t_plus) ? t_minus : t_plus;
			return true;
		}	

		return false;
	}
	
};

struct Plane{
	Vec p0, n; 
	Plane() : p0(Vec(0,0,0)), n(Vec(0,0,0)) {}
	Plane(const Vec & p0, const Vec& n) : p0(p0), n(n) {}

	Vec normal() const{ return p0-n; }
	bool intersect(const Ray &ray, double &t) const {
		Vec N = normal(); 
		double denom = dot(N, ray.d);
		if(denom > 1e-6){
			t = dot(p0-ray.e, N) / denom;
			return (t >= 0); 
		}
		return false;
	} 
};

struct Color{
	Vec color; 
	Color operator + (const Color& c) const { return Color(Vec(color.x+c.color.x, color.y+c.color.y, color.z+c.color.z)); }
	Color operator * (double d) const { return Color(Vec(color.x*d, color.y*d, color.z*d)); }
	Color() : color(Vec(0,0,0)){}
	Color(const Vec& rgb) : color(rgb){}
};

void col255(Color& col){
	col.color.x = (col.color.x > 255) ? 255 : (col.color.x < 0) ? 0 : col.color.x;
	col.color.y = (col.color.y > 255) ? 255 : (col.color.y < 0) ? 0 : col.color.y;
	col.color.z = (col.color.z > 255) ? 255 : (col.color.z < 0) ? 0 : col.color.z;
}


const int H = 500, W = 500;
const Sphere sphere1(Vec(W*0.5, H*0.5, 75), 50, 1); 
const Sphere sphere2(Vec(W*0.45, H*0.5, 25), 25, 0); 
const Sphere sphere3(Vec(W*0.5, H*0.9, 150), 10, 0); 
const Sphere sphere4(Vec(W*0.25, H*0.75, 100), 50, 0); 
const Sphere sphere5(Vec(W*0.75, H*0.75, 100), 50, 0); 
	const Color black(Vec(0, 0, 0)); //bkgd color
	const Color red(Vec(255, 0, 0));
	const Color green(Vec(0, 255, 0));
	const Color blue(Vec(0, 0, 255));
	const Color cyan(Vec(0, 255, 255)); //sphere1 color
	const Color white(Vec(255, 255, 255));
	const Color magenta(Vec(255, 0, 255)); 
	const Color yellow(Vec(255, 255, 0)); 
	const Color purple(Vec(102, 0, 102)); 	
void fillPixel(const Sphere &sphere, const Vec& light, Color& pixCol, Ray& ray, double& t, const Color& amb, const Color& dif, const Color& spec, double depth){
	if(sphere.intersect(ray,t)){
		if(depth < t) //check depth, don't fill pixel if something intersects ray first
			return;
		Vec p = ray.e + ray.d*t; //p(t) = e+td;		
		const Vec L = light - p;
		const Vec N = sphere.normal(p);
		const double ln = dot(L.normalize(), N.normalize());
		const Vec R = (N*ln)*2.0 - L; //reflection ray
		const Vec V = p - ray.e;
		const double vr = dot(V.normalize(), R.normalize()); 

		// I = I_a*k_a + I_d*k_d(l.n) + I_s*k_s(v.r)
		pixCol = ((amb+ dif*ln) * 0.5)+ ((spec*vr)*0.5);/// + (fillPixel(sphere, light, pixCol, R, t, amb, dif, spec, depth)*0.5); 
	/*	if(sphere.reflFlag == 1){
			std::cout << "hi got in" << std::endl;
			Ray reflRay = Ray(p, R); 
			double reflT;
			if(sphere1.intersect(reflRay, reflT)){
				pixCol = (pixCol + cyan)*0.5;
				std::cout << "reflected sphere1" << std::endl;
			}
			else if(sphere2.intersect(reflRay, reflT)){
				pixCol = (pixCol + white)*0.5;
				std::cout << "reflected sphere2" << std::endl;
			}
			else if(sphere3.intersect(reflRay, reflT)){
				pixCol = (pixCol + red)*0.5;
				std::cout << "reflected sphere3" << std::endl;
			}
			else if(sphere4.intersect(reflRay, reflT)){
				pixCol = (pixCol + green)*0.5; 
				std::cout << "reflected sphere4" << std::endl;
			}
			else if(sphere5.intersect(reflRay, reflT)){
				pixCol = (pixCol + blue)*0.5;
				std::cout << "reflected sphere5" << std::endl;
			}
		}*/
	}

}


void fillPlanePixel(const Plane &plane, const Vec& light, Color& pixCol, Ray& ray, double& t, const Color& amb, const Color& dif, const Color& spec, double depth){
	if(plane.intersect(ray,t)){
		if(depth < t) //check depth, don't fill pixel if something else intersects ray first
			return;
		Vec p = ray.e + ray.d*t; //p(t) = e+td;		
		//std::cout << p.x << ", " << p.y << ", " << p.z << std::endl;
		const Vec L = light - p;
		const Vec N = plane.normal();
		const double ln = dot(L.normalize(), N.normalize());
		const Vec R = (N*ln)*2.0 - L; //reflection ray
		const Vec V = p - ray.e;
		const double vr = dot(V.normalize(), R.normalize()); 


		// I = I_a*k_a + I_d*k_d(l.n) + I_s*k_s(v.r)
		pixCol = ((amb+ dif*ln) * 0.5)+ ((spec*vr)*0.5);/// + (fillPixel(sphere, light, pixCol, R, t, amb, dif, spec, depth)*0.5); 

		Ray shad_ray = Ray(p, L.normalize());    
		double shad_t; 


		//shadow cast
		if(sphere1.intersect(shad_ray, shad_t) || sphere2.intersect(shad_ray, shad_t) || sphere3.intersect(shad_ray, shad_t) || sphere4.intersect(shad_ray, shad_t) || sphere5.intersect(shad_ray, shad_t)){

			const Vec H = (L.normalize() + (ray.d*-1).normalize()).normalize();
			const double hn = dot(H, N.normalize()); 
			const Color black(Vec(0, 0, 0)); //bkgd color
			double max = 0; 
			if (ln > max)
				max = ln;
			pixCol = pixCol*0.5 + black*max*0.5 + black*hn*0.5;
			return;

		}


		
	}

}


const bool PROJ = true;

int main(){
	
	

	
	const Plane plane(Vec(0,500, 0), Vec(0,499,0)); 
	const Plane plane1(Vec(-1000, 250, 0), Vec(-999, 250, 0)); 
	const Plane plane2(Vec(1000, 250, 0), Vec(999, 250, 0)); 
//	const Plane plane3(Vec(250, 250, 1000), Vec(250, 250, 999)); 

	const Vec light(-100, -400, -200); //light point

	ofstream out("aalten.ppm");
	out << "P3\n";
	out << H << " " << W << "\n";
	out << "255\n"; 

	double t; //point on ray
	Color pixelCol(black); 
	
	for(int y = 0; y < H; y++){

		for(int x = 0; x < W; x++){				
			pixelCol = black;
			Color r(Vec(rand()%255, rand()%255, rand()%255));

			Ray ray;
			if (PROJ) {
				Vec c(W*0.5, H*0.5, -1000);
				Vec pt(x, y, 0);
				Vec unit_vector = (pt-c);
				ray = Ray(c, unit_vector);
	

			} else { //orthogonal
				ray = Ray(Vec(x,y,0), Vec(0,0,1));
			}
		
			double curr_t = 1000000;
		        fillPixel(sphere1, light, pixelCol, ray, t, cyan, cyan, white, curr_t);
			curr_t = t;
			fillPixel(sphere2, light, pixelCol, ray, t, r, r, white, curr_t);
			curr_t = t;
			fillPixel(sphere3, light, pixelCol, ray, t, red, red, white, curr_t);
			curr_t = t;
			fillPixel(sphere4, light, pixelCol, ray, t, green, green, white, curr_t);
			curr_t = t;
			fillPixel(sphere5, light, pixelCol, ray, t, blue, blue, white, curr_t);
			curr_t = t;
			fillPlanePixel(plane, light, pixelCol, ray, t, blue, blue, white, curr_t); 
			curr_t = t;
			/*fillPlanePixel(plane1, light, pixelCol, ray, t, magenta, magenta, white, curr_t); 
			curr_t = t;
			fillPlanePixel(plane2, light, pixelCol, ray, t, yellow, yellow, white, curr_t); 
			curr_t = t;*/
			//fillPlanePixel(plane3, light, pixelCol, ray, t, purple, purple, white, curr_t); 
			//curr_t = t;
			col255(pixelCol);
			
			
			out <<(int)pixelCol.color.x << ' ' << (int)pixelCol.color.y << ' ' << (int)pixelCol.color.z << '\n';
		}
	}



}
