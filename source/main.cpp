//////////////////////////////////////////////////////////////////////////////
//
//  --- main.cpp ---
//  Created by Brian Summa
//
//////////////////////////////////////////////////////////////////////////////

#include "common.h"
#include "SourcePath.h"
#include "common/lodepng.h"

using namespace Angel;

typedef vec4  color4;
typedef vec4  point4;


//Scene variables
enum { _SPHERE, _SQUARE, _BOX };
int scene = _SPHERE; //Simple sphere, square or cornell box
std::vector < Object* > sceneObjects;
point4 lightPosition;
color4 lightColor;
point4 cameraPosition;

//Recursion depth for raytracer
int maxDepth = 3;

namespace GLState {
    int window_width, window_height;

    bool render_line;

    std::vector < GLuint > objectVao;
    std::vector < GLuint > objectBuffer;

    GLuint vPosition, vNormal, vTexCoord;

    GLuint program;

    // Model-view and projection matrices uniform location
    GLuint  ModelView, ModelViewLight, NormalMatrix, Projection;

    //==========Trackball Variables==========
    static float curquat[4], lastquat[4];
    /* current transformation matrix */
    static float curmat[4][4];
    mat4 curmat_a;
    /* actual operation  */
    static int scaling;
    static int moving;
    static int panning;
    /* starting "moving" coordinates */
    static int beginx, beginy;
    /* ortho */
    float ortho_x, ortho_y;
    /* current scale factor */
    static float scalefactor;

    mat4  projection;
    mat4 sceneModelView;

    color4 light_ambient;
    color4 light_diffuse;
    color4 light_specular;

};


/* -------------------------------------------------------------------------- */
/* ----------------------  Write Image to Disk  ----------------------------- */
bool write_image(const char* filename, const unsigned char* Src,
    int Width, int Height, int channels) {
    unsigned bitdepth = 8;
    LodePNGColorType color_type;
    unsigned result;
    switch (channels)
    {
    case 1:
        color_type = LCT_GREY; break;
    case 2:
        color_type = LCT_GREY_ALPHA; break;
    case 3:
        color_type = LCT_RGB; break;
    default:
        color_type = LCT_RGBA; break;
    }
    result = lodepng_encode_file(filename, Src,
        static_cast<unsigned int>(Width), static_cast<unsigned int>(Height),
        color_type, bitdepth);
    if (result == 0)
        std::cerr << "finished writing " << filename << "." << std::endl;
    else
        std::cerr << "write to " << filename << " returned error code " << result << ". ("
        << lodepng_error_text(result) << ")" << std::endl;
    return result == 0;
}


/* -------------------------------------------------------------------------- */
/* -------- Given OpenGL matrices find ray in world coordinates of ---------- */
/* -------- window position x,y --------------------------------------------- */
std::vector < vec4 > findRay(GLdouble x, GLdouble y) {

    y = GLState::window_height - y;

    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    GLdouble modelViewMatrix[16];
    GLdouble projectionMatrix[16];
    for (unsigned int i = 0; i < 4; i++) {
        for (unsigned int j = 0; j < 4; j++) {
            modelViewMatrix[j * 4 + i] = GLState::sceneModelView[i][j];
            projectionMatrix[j * 4 + i] = GLState::projection[i][j];
        }
    }


    GLdouble nearPlaneLocation[3];
    _gluUnProject(x, y, 0.0, modelViewMatrix, projectionMatrix,
        viewport, &nearPlaneLocation[0], &nearPlaneLocation[1],
        &nearPlaneLocation[2]);

    GLdouble farPlaneLocation[3];
    _gluUnProject(x, y, 1.0, modelViewMatrix, projectionMatrix,
        viewport, &farPlaneLocation[0], &farPlaneLocation[1],
        &farPlaneLocation[2]);


    vec4 ray_origin = vec4(nearPlaneLocation[0], nearPlaneLocation[1], nearPlaneLocation[2], 1.0);
    vec3 temp = vec3(farPlaneLocation[0] - nearPlaneLocation[0],
        farPlaneLocation[1] - nearPlaneLocation[1],
        farPlaneLocation[2] - nearPlaneLocation[2]);
    temp = normalize(temp);
    vec4 ray_dir = vec4(temp.x, temp.y, temp.z, 0.0);

    std::vector < vec4 > result(2);
    result[0] = ray_origin;
    result[1] = ray_dir;

    return result;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
bool intersectionSort(Object::IntersectionValues i, Object::IntersectionValues j) {
    return (i.t_w < j.t_w);
}

/* -------------------------------------------------------------------------- */
/* ---------  Some debugging code: cast Ray = p0 + t*dir  ------------------- */
/* ---------  and print out what it hits =                ------------------- */
void castRayDebug(vec4 p0, vec4 dir) {

    std::vector < Object::IntersectionValues > intersections;

    for (unsigned int i = 0; i < sceneObjects.size(); i++) {
        intersections.push_back(sceneObjects[i]->intersect(p0, dir));
        intersections[intersections.size() - 1].ID_ = i;
    }

    for (unsigned int i = 0; i < intersections.size(); i++) {
        if (intersections[i].t_w != std::numeric_limits< double >::infinity()) {
            std::cout << "Hit " << intersections[i].name << " " << intersections[i].ID_ << "\n";
            std::cout << "P: " << intersections[i].P_w << "\n";
            std::cout << "N: " << intersections[i].N_w << "\n";
            vec4 L = lightPosition - intersections[i].P_w;
            L = normalize(L);
            std::cout << "L: " << L << "\n";
        }
    }

}


/* -------------------------------------------------------------------------- */
bool shadowFeeler(vec4 p0, Object* object) {
    bool inShadow = false;

    //TODO: Shadow code here
    float distance = sqrt((p0.x - lightPosition.x) * (p0.x - lightPosition.x) + (p0.y - lightPosition.y) * (p0.y - lightPosition.y) + (p0.z - lightPosition.z) * (p0.z - lightPosition.z));
    vec4 point_lightdir = normalize(lightPosition - p0);
    std::vector <float> tlist;
    for (int i = 0; i < sceneObjects.size(); i++) {
        tlist.push_back(sceneObjects[i]->intersect(p0, point_lightdir).t_w);
    }
    float tmin = tlist[0];
    int tminindex = 0;
    for (int i = 0; i < sceneObjects.size(); i++) {
        if (tlist[i] < tmin) {
            tmin = tlist[i];
            tminindex = i;
        }
    }
    if (tmin < distance) { inShadow = true; }
    return inShadow;
}

/* -------------------------------------------------------------------------- */
/* ----------  cast Ray = p0 + t*dir and intersect with sphere      --------- */
/* ----------  return color, right now shading is approx based      --------- */
/* ----------  depth                                                --------- */
bool CriticalAngleTest(float N12, vec4 N, vec4 I) {
    float critical = 1 - pow(N12, 2) * (1 - pow(dot(N, I), 2));
    if (critical > 0) {
        return false;
    }
    return true;
}
vec4 refract(vec4 dir, vec4 intersectpoint, vec4 N, float n1, float n2) {

    float n12 = n1 / n2;
    if (!CriticalAngleTest(n1 / n2, N, dir)) {
        return (-1 * n12 * dot(N, dir) - sqrt(1 - n12 * n12 * (1 - dot(N, dir) * dot(N, dir)))) * N + n12 * dir;
    }
    else { return vec4(0, 0, 0, 0); }
}

vec2 Fresnel(vec4 dir, vec4 intersectpoint, vec4 N, float n1, float n2) {
    float costheta_i = dot(dir, -N);
    float costheta_t = dot(refract(dir,intersectpoint,N,n1,n2),-N);
    float R_perp = pow((n1 * costheta_i - n2 * costheta_t) / (n1 * costheta_i + n2 * costheta_t),2);
    float R_para = pow((n2 * costheta_i - n1 * costheta_t) / (n2 * costheta_i + n1 * costheta_t),2);
    float R;
    
    if (!CriticalAngleTest((n1/n2),N,dir)) {
        R = (R_perp + R_para) / 2;
    }
    else {
        R = 1;
    }
    float T = 1 - R;

    return vec2(T, R);
}


vec4 castRay(vec4 p0, vec4 dir, Object* lastHitObject, int depth) {
    vec4 color = vec4(0.0, 0.0, 0.0, 0.0);

    if (depth > maxDepth) { return color; }

    //TODO: Raytracing code here
    std::vector < Object::IntersectionValues > intersections;

    for (unsigned int i = 0; i < sceneObjects.size(); i++) {
        intersections.push_back(sceneObjects[i]->intersect(p0, dir));
        intersections[intersections.size() - 1].ID_ = i;
    }
    int minindex = 0;
    float mintval = intersections[0].t_w;
    for (int i = 0; i < sceneObjects.size(); i++) {
        if (intersections[i].t_w < mintval) {
            mintval = intersections[i].t_w;
            minindex = i;
        }
    }
    if (mintval < EPSILON) {
        return color;
    }
    if (mintval == std::numeric_limits< double >::infinity()) {
        return color;
    }
    else {
        vec4 intersectpoint = intersections[minindex].P_w;
        if (sceneObjects[minindex]->shadingValues.Ks > 0) {
            vec4 R = reflect(dir, normalize(intersections[minindex].N_w));
            color += sceneObjects[minindex]->shadingValues.Ks * castRay(intersectpoint, R, lastHitObject, depth + 1);
        }
        
        if (sceneObjects[minindex]->shadingValues.Kt > 0) {
            vec4 R = reflect(dir, normalize(intersections[minindex].N_w));
            vec4 Normal = normalize(intersections[minindex].N_w);
            float n1 = 1;
            float n2 = 1.4;
            if (lastHitObject!=NULL) {
                n1 = 1.4;
                n2 = 1;
                Normal = -Normal;
                lastHitObject = NULL;
            }
            else { lastHitObject = sceneObjects[minindex]; }
            vec4 T = refract(dir, intersectpoint, Normal, n1, n2);
            vec2 TR = Fresnel(dir, intersectpoint, Normal, n1, n2);
            float KT = TR.x;
            float KR = TR.y;
            vec4 transmit_color = KT * castRay(intersectpoint, T, lastHitObject, depth + 1);
            vec4 reflect_color = KR * castRay(intersectpoint, R, lastHitObject, depth + 1);
            color += sceneObjects[minindex]->shadingValues.Kt * (transmit_color +reflect_color);
        }
        
        if (sceneObjects[minindex]->shadingValues.Kd > 0) {
            if (shadowFeeler(intersectpoint, lastHitObject)) {
                return vec4(0.0, 0.0, 0.0, 1.0);
            }
            else {
                //Compute phong shading
                vec4 L = normalize(lightPosition - intersectpoint);
                vec4 N = normalize(intersections[minindex].N_w);
                vec4 V = normalize(p0 - intersectpoint);
                vec4 R = normalize(reflect(-1 * L, N));

                float D = max(dot(L, N), 0.0);
                float S = pow(max(dot(V, R), 0.0), sceneObjects[minindex]->shadingValues.Kn);

                color += sceneObjects[minindex]->shadingValues.color * D * sceneObjects[minindex]->shadingValues.Kd;
                color += lightColor * S * sceneObjects[minindex]->shadingValues.Ks;
            }
        }
        color.x = fmin(color.x, 1.0);
        color.y = fmin(color.y, 1.0);
        color.z = fmin(color.z, 1.0);
        color.w = 1;
        return color;
    }
}



/* -------------------------------------------------------------------------- */
/* ------------  Ray trace our scene.  Output color to image and    --------- */
/* -----------   Output color to image and save to disk             --------- */
void rayTrace() {

    unsigned char* buffer = new unsigned char[GLState::window_width * GLState::window_height * 4];

    for (unsigned int i = 0; i < GLState::window_width; i++) {
        for (unsigned int j = 0; j < GLState::window_height; j++) {

            int idx = j * GLState::window_width + i;
            std::vector < vec4 > ray_o_dir = findRay(i, j);
            vec4 color = castRay(ray_o_dir[0], vec4(ray_o_dir[1].x, ray_o_dir[1].y, ray_o_dir[1].z, 0.0), NULL, 0);
            buffer[4 * idx] = color.x * 255;
            buffer[4 * idx + 1] = color.y * 255;
            buffer[4 * idx + 2] = color.z * 255;
            buffer[4 * idx + 3] = color.w * 255;
        }
    }

    write_image("output.png", buffer, GLState::window_width, GLState::window_height, 4);

    delete[] buffer;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
static void error_callback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void initGL() {

    GLState::light_ambient = vec4(lightColor.x, lightColor.y, lightColor.z, 1.0);
    GLState::light_diffuse = vec4(lightColor.x, lightColor.y, lightColor.z, 1.0);
    GLState::light_specular = vec4(lightColor.x, lightColor.y, lightColor.z, 1.0);


    std::string vshader = source_path + "/shaders/vshader.glsl";
    std::string fshader = source_path + "/shaders/fshader.glsl";

    GLchar* vertex_shader_source = readShaderSource(vshader.c_str());
    GLchar* fragment_shader_source = readShaderSource(fshader.c_str());

    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, (const GLchar**)&vertex_shader_source, NULL);
    glCompileShader(vertex_shader);
    check_shader_compilation(vshader, vertex_shader);

    GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, (const GLchar**)&fragment_shader_source, NULL);
    glCompileShader(fragment_shader);
    check_shader_compilation(fshader, fragment_shader);

    GLState::program = glCreateProgram();
    glAttachShader(GLState::program, vertex_shader);
    glAttachShader(GLState::program, fragment_shader);

    glLinkProgram(GLState::program);
    check_program_link(GLState::program);

    glUseProgram(GLState::program);

    glBindFragDataLocation(GLState::program, 0, "fragColor");

    // set up vertex arrays
    GLState::vPosition = glGetAttribLocation(GLState::program, "vPosition");
    GLState::vNormal = glGetAttribLocation(GLState::program, "vNormal");

    // Retrieve transformation uniform variable locations
    GLState::ModelView = glGetUniformLocation(GLState::program, "ModelView");
    GLState::NormalMatrix = glGetUniformLocation(GLState::program, "NormalMatrix");
    GLState::ModelViewLight = glGetUniformLocation(GLState::program, "ModelViewLight");
    GLState::Projection = glGetUniformLocation(GLState::program, "Projection");

    if (GLState::objectVao.size() > 0) {
        glDeleteVertexArrays(GLState::objectVao.size(), &GLState::objectVao[0]);
    }
    GLState::objectVao.resize(sceneObjects.size());
    glGenVertexArrays(sceneObjects.size(), &GLState::objectVao[0]);

    if (GLState::objectBuffer.size() > 0) {
        glDeleteBuffers(GLState::objectBuffer.size(), &GLState::objectBuffer[0]);
    }
    GLState::objectBuffer.resize(sceneObjects.size());
    glGenBuffers(sceneObjects.size(), &GLState::objectBuffer[0]);

    for (unsigned int i = 0; i < sceneObjects.size(); i++) {
        glBindVertexArray(GLState::objectVao[i]);
        glBindBuffer(GL_ARRAY_BUFFER, GLState::objectBuffer[i]);
        size_t vertices_bytes = sceneObjects[i]->mesh.vertices.size() * sizeof(vec4);
        size_t normals_bytes = sceneObjects[i]->mesh.normals.size() * sizeof(vec3);

        glBufferData(GL_ARRAY_BUFFER, vertices_bytes + normals_bytes, NULL, GL_STATIC_DRAW);
        size_t offset = 0;
        glBufferSubData(GL_ARRAY_BUFFER, offset, vertices_bytes, &sceneObjects[i]->mesh.vertices[0]);
        offset += vertices_bytes;
        glBufferSubData(GL_ARRAY_BUFFER, offset, normals_bytes, &sceneObjects[i]->mesh.normals[0]);

        glEnableVertexAttribArray(GLState::vNormal);
        glEnableVertexAttribArray(GLState::vPosition);

        glVertexAttribPointer(GLState::vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0));
        glVertexAttribPointer(GLState::vNormal, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(vertices_bytes));

    }



    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);

    glClearColor(0.8, 0.8, 1.0, 1.0);

    //Quaternion trackball variables, you can ignore
    GLState::scaling = 0;
    GLState::moving = 0;
    GLState::panning = 0;
    GLState::beginx = 0;
    GLState::beginy = 0;

    TrackBall::matident(GLState::curmat);
    TrackBall::trackball(GLState::curquat, 0.0f, 0.0f, 0.0f, 0.0f);
    TrackBall::trackball(GLState::lastquat, 0.0f, 0.0f, 0.0f, 0.0f);
    TrackBall::add_quats(GLState::lastquat, GLState::curquat, GLState::curquat);
    TrackBall::build_rotmatrix(GLState::curmat, GLState::curquat);

    GLState::scalefactor = 1.0;
    GLState::render_line = false;

}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void initCornellBox() {
    cameraPosition = point4(0.0, 0.0, 6.0, 1.0);
    lightPosition = point4(0.0, 1.5, 0.0, 1.0);
    lightColor = color4(1.0, 1.0, 1.0, 1.0);

    sceneObjects.clear();

    { //Back Wall
        sceneObjects.push_back(new Square("Back Wall"));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(1.0, 1.0, 1.0, 1.0);
        _shadingValues.Ka = 0.0;
        _shadingValues.Kd = 1.0;
        _shadingValues.Ks = 0.0;
        _shadingValues.Kn = 16.0;
        _shadingValues.Kt = 0.0;
        _shadingValues.Kr = 0.0;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(Translate(0.0, 0.0, -2.0) * Scale(2.0, 2.0, 1.0));
    }

    { //Left Wall
        sceneObjects.push_back(new Square("Left Wall"));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(1.0, 0.0, 0.0, 1.0);
        _shadingValues.Ka = 0.0;
        _shadingValues.Kd = 1.0;
        _shadingValues.Ks = 0.0;
        _shadingValues.Kn = 16.0;
        _shadingValues.Kt = 0.0;
        _shadingValues.Kr = 0.0;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(RotateY(90) * Translate(0.0, 0.0, -2.0) * Scale(2.0, 2.0, 1.0));
    }

    { //Right Wall
        sceneObjects.push_back(new Square("Right Wall"));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(0.5, 0.0, 0.5, 1.0);
        _shadingValues.Ka = 0.0;
        _shadingValues.Kd = 1.0;
        _shadingValues.Ks = 0.0;
        _shadingValues.Kn = 16.0;
        _shadingValues.Kt = 0.0;
        _shadingValues.Kr = 0.0;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(RotateY(-90) * Translate(0.0, 0.0, -2.0) * Scale(2.0, 2.0, 1.0));
    }

    { //Floor
        sceneObjects.push_back(new Square("Floor"));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(1.0, 1.0, 1.0, 1.0);
        _shadingValues.Ka = 0.0;
        _shadingValues.Kd = 1.0;
        _shadingValues.Ks = 0.0;
        _shadingValues.Kn = 16.0;
        _shadingValues.Kt = 0.0;
        _shadingValues.Kr = 0.0;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(RotateX(-90) * Translate(0.0, 0.0, -2.0) * Scale(2.0, 2.0, 1.0));
    }

    { //Ceiling
        sceneObjects.push_back(new Square("Ceiling"));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(1.0, 1.0, 1.0, 1.0);
        _shadingValues.Ka = 0.0;
        _shadingValues.Kd = 1.0;
        _shadingValues.Ks = 0.0;
        _shadingValues.Kn = 16.0;
        _shadingValues.Kt = 0.0;
        _shadingValues.Kr = 0.0;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(RotateX(90) * Translate(0.0, 0.0, -2.0) * Scale(2.0, 2.0, 1.0));
    }

    { //Front Wall
        sceneObjects.push_back(new Square("Front Wall"));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(1.0, 1.0, 1.0, 1.0);
        _shadingValues.Ka = 0.0;
        _shadingValues.Kd = 1.0;
        _shadingValues.Ks = 0.0;
        _shadingValues.Kn = 16.0;
        _shadingValues.Kt = 0.0;
        _shadingValues.Kr = 0.0;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(RotateY(180) * Translate(0.0, 0.0, -2.0) * Scale(2.0, 2.0, 1.0));
    }


    {
        sceneObjects.push_back(new Sphere("Glass sphere"));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(1.0, 0.0, 0.0, 1.0);
        _shadingValues.Ka = 0.0;
        _shadingValues.Kd = 0.0;
        _shadingValues.Ks = 0.0;
        _shadingValues.Kn = 16.0;
        _shadingValues.Kt = 1.0;
        _shadingValues.Kr = 1.4;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(Translate(1.0, -1.25, 0.5) * Scale(0.75, 0.75, 0.75));
    }

    {
        sceneObjects.push_back(new Sphere("Mirrored Sphere"));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(1.0, 1.0, 1.0, 1.0);
        _shadingValues.Ka = 0.0;
        _shadingValues.Kd = 0.0;
        _shadingValues.Ks = 1.0;
        _shadingValues.Kn = 16.0;
        _shadingValues.Kt = 0.0;
        _shadingValues.Kr = 0.0;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(Translate(-1.0, -1.25, -1.0) * Scale(0.75, 0.75, 0.75));
    }

    initGL();
}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void initUnitSphere() {
    cameraPosition = point4(0.0, 0.0, 3.0, 1.0);
    lightPosition = point4(0.0, 0.0, 4.0, 1.0);
    lightColor = color4(1.0, 1.0, 1.0, 1.0);

    sceneObjects.clear();

    {
        sceneObjects.push_back(new Sphere("Diffuse sphere"));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(1.0, 0.0, 0.0, 1.0);
        _shadingValues.Ka = 0.0;
        _shadingValues.Kd = 1.0;
        _shadingValues.Ks = 0.0;
        _shadingValues.Kn = 16.0;
        _shadingValues.Kt = 0.0;
        _shadingValues.Kr = 0.0;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(Translate(1, 0, 0));

        sceneObjects.push_back(new Sphere("Diffuse sphere"));
        Object::ShadingValues _shadingValues2;
        _shadingValues2.color = vec4(0.0, 1.0, 0.0, 1.0);
        _shadingValues2.Ka = 0.0;
        _shadingValues2.Kd = 1.0;
        _shadingValues2.Ks = 0.0;
        _shadingValues2.Kn = 16.0;
        _shadingValues2.Kt = 0.0;
        _shadingValues2.Kr = 0.0;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues2);
        sceneObjects[sceneObjects.size() - 1]->setModelView(Translate(-1, 0, 0) * Scale(1.5, 1.5, 1.5));
    }

    initGL();
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void initUnitSquare() {
    cameraPosition = point4(0.0, 0.0, 3.0, 1.0);
    lightPosition = point4(0.0, 0.0, 4.0, 1.0);
    lightColor = color4(1.0, 1.0, 1.0, 1.0);

    sceneObjects.clear();

    { //Back Wall
        sceneObjects.push_back(new Square("Unit Square"));
        Object::ShadingValues _shadingValues;
        _shadingValues.color = vec4(1.0, 0.0, 0.0, 1.0);
        _shadingValues.Ka = 0.0;
        _shadingValues.Kd = 1.0;
        _shadingValues.Ks = 0.0;
        _shadingValues.Kn = 16.0;
        _shadingValues.Kt = 0.0;
        _shadingValues.Kr = 0.0;
        sceneObjects[sceneObjects.size() - 1]->setShadingValues(_shadingValues);
        sceneObjects[sceneObjects.size() - 1]->setModelView(mat4());
    }

    initGL();

}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    if (key == GLFW_KEY_1 && action == GLFW_PRESS) {
        scene = _SPHERE;
        initUnitSphere();
    }
    if (key == GLFW_KEY_2 && action == GLFW_PRESS) {
        scene = _SQUARE;
        initUnitSquare();
    }
    if (key == GLFW_KEY_3 && action == GLFW_PRESS) {
        scene = _BOX;
        initCornellBox();
    }
    if (key == GLFW_KEY_R && action == GLFW_PRESS)
        rayTrace();
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
static void mouseClick(GLFWwindow* window, int button, int action, int mods) {

    if (GLFW_RELEASE == action) {
        GLState::moving = GLState::scaling = GLState::panning = false;
        return;
    }

    if (mods & GLFW_MOD_SHIFT) {
        GLState::scaling = true;
    }
    else if (mods & GLFW_MOD_ALT) {
        GLState::panning = true;
    }
    else {
        GLState::moving = true;
        TrackBall::trackball(GLState::lastquat, 0, 0, 0, 0);
    }

    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    GLState::beginx = xpos; GLState::beginy = ypos;

    std::vector < vec4 > ray_o_dir = findRay(xpos, ypos);
    castRayDebug(ray_o_dir[0], vec4(ray_o_dir[1].x, ray_o_dir[1].y, ray_o_dir[1].z, 0.0));

}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void mouseMove(GLFWwindow* window, double x, double y) {

    int W, H;
    glfwGetFramebufferSize(window, &W, &H);


    float dx = (x - GLState::beginx) / (float)W;
    float dy = (GLState::beginy - y) / (float)H;

    if (GLState::panning)
    {
        GLState::ortho_x += dx;
        GLState::ortho_y += dy;

        GLState::beginx = x; GLState::beginy = y;
        return;
    }
    else if (GLState::scaling)
    {
        GLState::scalefactor *= (1.0f + dx);

        GLState::beginx = x; GLState::beginy = y;
        return;
    }
    else if (GLState::moving)
    {
        TrackBall::trackball(GLState::lastquat,
            (2.0f * GLState::beginx - W) / W,
            (H - 2.0f * GLState::beginy) / H,
            (2.0f * x - W) / W,
            (H - 2.0f * y) / H
        );

        TrackBall::add_quats(GLState::lastquat, GLState::curquat, GLState::curquat);
        TrackBall::build_rotmatrix(GLState::curmat, GLState::curquat);

        GLState::beginx = x; GLState::beginy = y;
        return;
    }
}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void drawObject(Object* object, GLuint vao, GLuint buffer) {

    color4 material_ambient(object->shadingValues.color.x * object->shadingValues.Ka,
        object->shadingValues.color.y * object->shadingValues.Ka,
        object->shadingValues.color.z * object->shadingValues.Ka, 1.0);
    color4 material_diffuse(object->shadingValues.color.x,
        object->shadingValues.color.y,
        object->shadingValues.color.z, 1.0);
    color4 material_specular(object->shadingValues.Ks,
        object->shadingValues.Ks,
        object->shadingValues.Ks, 1.0);
    float  material_shininess = object->shadingValues.Kn;

    color4 ambient_product = GLState::light_ambient * material_ambient;
    color4 diffuse_product = GLState::light_diffuse * material_diffuse;
    color4 specular_product = GLState::light_specular * material_specular;

    glUniform4fv(glGetUniformLocation(GLState::program, "AmbientProduct"), 1, ambient_product);
    glUniform4fv(glGetUniformLocation(GLState::program, "DiffuseProduct"), 1, diffuse_product);
    glUniform4fv(glGetUniformLocation(GLState::program, "SpecularProduct"), 1, specular_product);
    glUniform4fv(glGetUniformLocation(GLState::program, "LightPosition"), 1, lightPosition);
    glUniform1f(glGetUniformLocation(GLState::program, "Shininess"), material_shininess);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
    glVertexAttribPointer(GLState::vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0));
    glVertexAttribPointer(GLState::vNormal, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(object->mesh.vertices.size() * sizeof(vec4)));

    mat4 objectModelView = GLState::sceneModelView * object->getModelView();


    glUniformMatrix4fv(GLState::ModelViewLight, 1, GL_TRUE, GLState::sceneModelView);
    glUniformMatrix3fv(GLState::NormalMatrix, 1, GL_TRUE, Normal(objectModelView));
    glUniformMatrix4fv(GLState::ModelView, 1, GL_TRUE, objectModelView);

    glDrawArrays(GL_TRIANGLES, 0, object->mesh.vertices.size());

}


int main(void) {

    GLFWwindow* window;

    glfwSetErrorCallback(error_callback);

    if (!glfwInit())
        exit(EXIT_FAILURE);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    glfwWindowHint(GLFW_SAMPLES, 4);

    window = glfwCreateWindow(768, 768, "Raytracer", NULL, NULL);
    if (!window) {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwSetKeyCallback(window, keyCallback);
    glfwSetMouseButtonCallback(window, mouseClick);
    glfwSetCursorPosCallback(window, mouseMove);


    glfwMakeContextCurrent(window);
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
    glfwSwapInterval(1);

    switch (scene) {
    case _SPHERE:
        initUnitSphere();
        break;
    case _SQUARE:
        initUnitSquare();
        break;
    case _BOX:
        initCornellBox();
        break;
    }

    while (!glfwWindowShouldClose(window)) {

        int width, height;
        glfwGetFramebufferSize(window, &width, &height);

        GLState::window_height = height;
        GLState::window_width = width;

        glViewport(0, 0, width, height);


        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        mat4 track_ball = mat4(GLState::curmat[0][0], GLState::curmat[1][0],
            GLState::curmat[2][0], GLState::curmat[3][0],
            GLState::curmat[0][1], GLState::curmat[1][1],
            GLState::curmat[2][1], GLState::curmat[3][1],
            GLState::curmat[0][2], GLState::curmat[1][2],
            GLState::curmat[2][2], GLState::curmat[3][2],
            GLState::curmat[0][3], GLState::curmat[1][3],
            GLState::curmat[2][3], GLState::curmat[3][3]);

        GLState::sceneModelView = Translate(-cameraPosition) *   //Move Camera Back
            Translate(GLState::ortho_x, GLState::ortho_y, 0.0) *
            track_ball *                   //Rotate Camera
            Scale(GLState::scalefactor,
                GLState::scalefactor,
                GLState::scalefactor);   //User Scale

        GLfloat aspect = GLfloat(width) / height;

        switch (scene) {
        case _SPHERE:
        case _SQUARE:
            GLState::projection = Perspective(45.0, aspect, 0.01, 100.0);
            break;
        case _BOX:
            GLState::projection = Perspective(45.0, aspect, 4.5, 100.0);
            break;
        }

        glUniformMatrix4fv(GLState::Projection, 1, GL_TRUE, GLState::projection);

        for (unsigned int i = 0; i < sceneObjects.size(); i++) {
            drawObject(sceneObjects[i], GLState::objectVao[i], GLState::objectBuffer[i]);
        }

        glfwSwapBuffers(window);
        glfwPollEvents();

    }

    glfwDestroyWindow(window);

    glfwTerminate();
    exit(EXIT_SUCCESS);
}
