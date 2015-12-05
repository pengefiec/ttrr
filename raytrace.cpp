//
// template-rt.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
int g_width;
int g_height;

struct Ray
{
    vec4 origin;
    vec4 dir;
};
// TODO: add structs for spheres, lights and anything else you may need.
struct Sphere
{
    string name;
    vec4 pos;
    float scale_x;
    float scale_y;
    float scale_z;
    float red;
    float green;
    float blue;
    float ka;
    float kd;
    float ks;
    float kr;
    float n;

};
struct Light
{
    string name;
    vec4 pos;
    float ir;
    float ig;
    float ib;
};
struct Background
{
    float red;
    float green;
    float blue;
};
struct Ambient
{
    float ir;
    float ig;
    float ib;
};
struct Intersection_key
{
    vec4 origin;
    Ray ray;
};
string output;
/*all colors of pixel*/
vector<vec4> g_colors;
/*all intersections*/
vector<vec4> intersections;
/*all spheres*/
vector<Sphere> spheres;
/*all lights*/
vector<Light> lights;

Ambient ambient;
Background back;

//what these used for?
float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;


// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    //TODO: add parsing of NEAR, LEFT, RIGHT, BOTTOM, TOP, SPHERE, LIGHT, BACK, AMBIENT, OUTPUT.
    const int num_lables = 11;
    const string labels[] = {"NEAR", "LEFT", "RIGHT", "BOTTOM", "TOP", "RES", "SPHERE", "LIGHT", "BACK", "AMBIENT", "OUTPUT"};
    unsigned label_id = find(labels, labels + num_lables, vs[0]) - labels;
    switch (label_id)
    {
        case 0: g_near = toFloat(vs[1]);
            break;
        case 1: g_left = toFloat(vs[1]);
            break;
        case 2: g_right = toFloat(vs[1]);
            break;
        case 3: g_bottom = toFloat(vs[1]);
            break;
        case 4: g_top = toFloat(vs[1]);
            break;
        case 5: g_width = (int)toFloat(vs[1]);
                g_height = (int)toFloat(vs[2]);
                g_colors.resize(g_width * g_height);
            break;
        case 6: {
            Sphere s;
            s.name = vs[1];
            s.pos = vec4(toFloat(vs[2]), toFloat(vs[3]), toFloat(vs[4]), 1.0);
            s.scale_x = toFloat(vs[5]);
            s.scale_y = toFloat(vs[6]);
            s.scale_z = toFloat(vs[7]);
            s.red = toFloat(vs[8]);
            s.green = toFloat(vs[9]);
            s.blue = toFloat(vs[10]);
            s.ka = toFloat(vs[11]);
            s.kd = toFloat(vs[12]);
            s.ks = toFloat(vs[13]);
            //what kr used for
            s.kr = toFloat(vs[14]);
            s.n = toFloat(vs[15]);
            spheres.push_back(s);
            break;
        }
        case 7: {
            Light l;
            l.name = vs[1];
            l.pos = vec4(toFloat(vs[2]), toFloat(vs[3]), toFloat(vs[4]), 1.0);
            l.ir = toFloat(vs[5]);
            l.ig = toFloat(vs[6]);
            l.ib = toFloat(vs[7]);
            lights.push_back(l);
            break;
        }
        case 8: {
            back.red = toFloat(vs[1]);
            back.green = toFloat(vs[2]);
            back.blue = toFloat(vs[3]);
            break;
        }
        case 9: {
            ambient.ir = toFloat(vs[1]);
            ambient.ig = toFloat(vs[2]);
            ambient.ib = toFloat(vs[3]);
            break;
        }
        case 10:
            output = vs[1];
    }
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
    //cout<<"ix is " << ix << "iy is " << iy << endl;
    //printv(color);
}



/*  Intersection routine
 *  the start_point could be the intersection point on a sphere surface
 *  return the delta value, if delta >= 0 and in the view volum, then it's valid intersection
 *  return the delta value, if delta >= 0 and in the view volum, then it's valid intersection
 */
// TODO: add your ray-sphere intersection routine here.
 bool intersect(const Sphere &sphere, const Ray &ray, const vec4 &start_point, vec4 &intersect_point, float min_dist, bool init) {
    mat4 transform = mat4();
    mat4 translate = Translate(sphere.pos[0], sphere.pos[1], sphere.pos[2]);
    mat4 scale = Scale(sphere.scale_x, sphere.scale_y, sphere.scale_z);
    transform = transform * translate * scale;
    mat4 i_transform;
    InvertMatrix(transform, i_transform);
    vec4 new_origin = i_transform * start_point;
    vec4 new_dir = i_transform * ray.dir;
    float sx = new_origin[0];
    float sy = new_origin[1];
    float sz = new_origin[2];
    float cx = new_dir[0];
    float cy = new_dir[1];
    float cz = new_dir[2];
    float a = cx * cx + cy * cy + cz * cz;
    float b = 2 * (sx * cx + sy * cy + sz * cz);
    float c = sx * sx + sy * sy + sz * sz - 1;
    float delta = b * b - 4 * a * c;
    if (delta >= 0) {
        float hit_1 = (-1 * b + sqrt(delta)) / (2 * a);
        float hit_2 = (-1 * b - sqrt(delta)) / (2 * a);
        float hit;
        float hit_1_dist = abs(start_point.z + ray.dir.z * hit_1);
        float hit_2_dist = abs(start_point.z + ray.dir.z * hit_2);
        if(hit_1 < 0 && hit_2 < 0 && init)
            return false;
//        else if(hit_1 < 0)
//            hit = hit_2;
//        else if(hit_2 < 0)
//            hit = hit_1;
        else if ( hit_1_dist < min_dist && hit_2_dist < min_dist)
            return false;
        else if(hit_1_dist < min_dist && hit_2_dist > min_dist)
            hit = hit_2;
        else if(hit_1_dist > min_dist && hit_2_dist < min_dist)
            hit = hit_1;
        else if(hit_1 > hit_2)
            hit = hit_2;
        else
            hit = hit_1;

        if ( hit_1 < min_dist && hit_2 < min_dist)
            return false;
        else if(hit_1 < min_dist && hit_2 > min_dist)
            hit = hit_2;
        else if(hit_1 > min_dist && hit_2 < min_dist)
            hit = hit_1;
        else if(hit_1 > hit_2)
            hit = hit_2;
        else
            hit = hit_1;
        intersect_point = vec4(start_point[0] + ray.dir[0] * hit, start_point[1] + ray.dir[1] * hit, start_point[2] + ray.dir[2] * hit, 1.0f);
        //printv(intersect_point);
        return true;
    }
    return false;
}

/*
 * calculate the distance between a point and the ray origin
 */
float calDist(const Ray &ray, const vec4 &point) {
    float dx = point[0] - ray.origin[0];
    float dy = point[1] - ray.origin[1];
    float dz = point[2] - ray.origin[2];
    return sqrt(dx * dx + dy * dy + dz * dz);
}

/*
 * trace light
 */
vec3 lightTrace(const vec4 &point, const Light &light, const Sphere &sphere, const vec4 &view_dir) {
    // the normal of a point on sphere surface
    //vec4 normal = normalize(point - sphere.pos);
    vec4 rawNormal =  point - sphere.pos;
    vec4 normal = normalize(vec4(rawNormal.x / (sphere.scale_x * sphere.scale_x), rawNormal.y / (sphere.scale_y * sphere.scale_y), rawNormal.z / (sphere.scale_z * sphere.scale_z), 0.0f));
    //to calculate the theta, the light direction need to point to the light
    vec4 light_dir = normalize(light.pos - point);
    //the angle between normal and light direction
    float cos_norm_light = dot(normal, light_dir);
    Ray r_ray;
    r_ray.origin = point;
    r_ray.dir = light_dir;
    vec4 intersect_point = vec4();
    bool isIntersect = false;
    for(unsigned i = 0; i < spheres.size(); ++i){
        if(spheres[i].name.compare(sphere.name) != 0) {
            if(intersect(spheres[i], r_ray, point, intersect_point, 0.001f, false))
                isIntersect = true;
        }
    }

    float ir = 0, ig = 0, ib = 0;
    if(cos_norm_light > 0 && !isIntersect) {
        //get the reflect light direction
        vec4 reflect_dir = normalize(light_dir - 2 * cos_norm_light * normal);
        //compute the reflection color for r, g, b separately
        vec4 vd = normalize(-1 * view_dir);
        float rdotv = dot(vd, reflect_dir);
        if(rdotv > 0){
            rdotv = 0;
        }
        float cos_view_reflect_n = pow(rdotv, sphere.n);
        ir = light.ir * sphere.kd * cos_norm_light * sphere.red + light.ir * sphere.ks * cos_view_reflect_n;
        ig = light.ig * sphere.kd * cos_norm_light * sphere.green + light.ig * sphere.ks * cos_view_reflect_n;
        ib = light.ib * sphere.kd * cos_norm_light * sphere.blue + light.ib * sphere.ks * cos_view_reflect_n;
    }
    return vec3(ir, ig, ib);
}

/*
 * calculate the lighting
 */
vec3 calLighting(const vec4 &point, const Sphere &sphere, const Ray &ray) {
    vec3 light_color_sum = vec3();
    for(unsigned i = 0; i < lights.size(); ++i) {
//        Ray ray;
//        ray.origin = point;
//        ray.dir = normalize(lights[i].pos - point);
//        printv(ray.dir);
//        bool isIntersect = intersect(sphere, ray, point, intersect_point, 0.0001);
       //if(!isIntersect){
        light_color_sum += lightTrace(point, lights[i], sphere, ray.dir);
       //}
    }
    return light_color_sum;
}

/*
 * Ray tracing
 */
vec3 rayTrace(const Ray& ray, int counter)
{
    if(--counter == -1) {
        return vec3();
    }
    // TODO: implement your ray tracing routine here.
    float min_dist = MAXFLOAT;
    bool isIntersect;
    vec4 closest_point = vec4();
    int intersect_sphere_index;

    for(unsigned i = 0; i < spheres.size(); ++i) {
        vec4 intersect_point = vec4();
        float threshold = counter == 2 ? 1.0f : 0.001f;
        if(intersect(spheres[i],ray, ray.origin, intersect_point, threshold, counter == 2)){
            float dist = calDist(ray, intersect_point);
            if(dist < min_dist) {
                min_dist = dist;
                closest_point = intersect_point;
                intersect_sphere_index = i;
                isIntersect = true;
            }

        }
    }
    if(counter < 2 && !isIntersect){
        return vec3();
    }
    if(!isIntersect) {
        return vec3(back.red, back.green, back.blue);
    }

    Sphere closest_sphere = spheres[intersect_sphere_index];
    vec3 color_local = calLighting(closest_point, closest_sphere, ray);
    vec3 color_ambient = closest_sphere.ka * vec3(ambient.ir * closest_sphere.red, ambient.ig * closest_sphere.green,
                                                  ambient.ib * closest_sphere.blue);
    vec4 normal = normalize(closest_point - spheres[intersect_sphere_index].pos);
    vec4 reflect_dir = normalize(ray.dir - 2 * (dot(ray.dir, normal) * normal ));
    Ray new_ray;
    new_ray.origin = closest_point;
    new_ray.dir = reflect_dir;
    vec3 color_reflect = closest_sphere.kr * rayTrace(new_ray, counter);
    //printv(color_reflect);
    return color_ambient + color_local + color_reflect;
}
vec4 trace(const Ray &ray) {
    vec3 color = rayTrace(ray, 3);
    return vec4(color[0], color[1], color[2], 1.0f);
}
/*
 * get the direction of the ray from pixel ix, iy of the near plane
 */

vec4 getDir(int ix, int iy)
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
    float h = (g_top - g_bottom)/2;
    float w = (g_right - g_left)/2;
    //calculate the pixel coordinate in camera coordinate system
    float ux = -1 * h + h * 2 * ix / (g_height - 1);
    float uy = -1 * w + w * 2 * iy / (g_width - 1);
    vec4 dir;
    dir = normalize (vec4(ux, uy, -1 * g_near, 0.0f));
    return dir;
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray);
    //printv(color);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, char* fname, unsigned char* pixels) 
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }
    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // TODO: clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++){
                float color = ((float*)g_colors[y*g_width+x])[i] > 1.0f ? 1.0f : ((float*)g_colors[y*g_width+x])[i];
                buf[y*g_width*3+x*3+i] = (unsigned char)(color * 255.9f);
            }

    // TODO: change file name based on input file name.
    savePPM(g_width, g_height, &output[0], buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    //string param = "../assignment3testsandresults/testImgPlane.txt";
    //char *paramc = &param[0];
    loadFile(argv[1]);
    render();
    saveFile();
	return 0;
}

