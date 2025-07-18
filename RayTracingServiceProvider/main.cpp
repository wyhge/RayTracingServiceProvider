#include <SDL3\SDL.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <random>

const int WIDTH = 800;
const int HEIGHT = 600;
const int MAX_DEPTH = 5;
const float PI = 3.14159265358979323846f;

// 向量类
struct Vec3 {
    float x, y, z;

    Vec3(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}

    Vec3 operator + (const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator - (const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator * (float f) const { return Vec3(x * f, y * f, z * f); }
    Vec3 operator / (float f) const { return Vec3(x / f, y / f, z / f); }
    Vec3 mult(const Vec3& v) const { return Vec3(x * v.x, y * v.y, z * v.z); } // 分量乘法

    float length() const { return std::sqrt(x * x + y * y + z * z); }
    Vec3 normalize() const { return (*this) / length(); }
    float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3 cross(const Vec3& v) const {
        return Vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
};

// 光线类
struct Ray {
    Vec3 origin;
    Vec3 direction;

    Ray(const Vec3& o, const Vec3& d) : origin(o), direction(d) {}

    Vec3 pointAt(float t) const {
        return origin + direction * t;
    }
};

// 材质类型
enum MaterialType { DIFFUSE, METAL, GLASS };

// 材质类
struct Material {
    MaterialType type;
    Vec3 albedo;        // 颜色
    float roughness;    // 金属粗糙度 (0-1)
    float ior;          // 折射率 (玻璃用)

    Material(MaterialType t, const Vec3& a, float r = 0.0f, float i = 1.5f)
        : type(t), albedo(a), roughness(r), ior(i) {
    }
};

// 球体类
struct Sphere {
    Vec3 center;
    float radius;
    Material material;

    Sphere(const Vec3& c, float r, const Material& m)
        : center(c), radius(r), material(m) {
    }

    // 光线与球体求交
    bool intersect(const Ray& ray, float& t) const {
        Vec3 oc = ray.origin - center;
        float a = ray.direction.dot(ray.direction);
        float b = 2.0f * oc.dot(ray.direction);
        float c = oc.dot(oc) - radius * radius;
        float discriminant = b * b - 4 * a * c;

        if (discriminant < 0) return false;

        float sqrtd = std::sqrt(discriminant);
        float t0 = (-b - sqrtd) / (2.0f * a);
        float t1 = (-b + sqrtd) / (2.0f * a);

        if (t0 > 0) {
            t = t0;
            return true;
        }
        else if (t1 > 0) {
            t = t1;
            return true;
        }
        return false;
    }
};

// 相机类
struct Camera {
    Vec3 position;      // 相机位置
    Vec3 forward;       // 前方向
    Vec3 up;            // 上方向
    Vec3 right;         // 右方向
    float fov;          // 视野角度 (度)

    Camera(const Vec3& pos, const Vec3& target, float fovDeg)
        : position(pos), fov(fovDeg)
    {
        forward = (target - position).normalize();
        right = forward.cross(Vec3(0, 1, 0)).normalize(); // 假设世界Y轴为上
        up = right.cross(forward).normalize();
    }

    // 生成光线
    Ray getRay(float u, float v) const {
        // 转换为弧度
        float aspect = float(WIDTH) / HEIGHT;
        float fovScale = std::tan(fov * PI / 360.0f);

        // 计算屏幕上的点
        Vec3 screenPoint =
            right * (u * aspect * fovScale) +
            up * (v * fovScale) +
            position + forward;

        return Ray(position, (screenPoint - position).normalize());
    }
};

// 场景类
struct Scene {
    std::vector<Sphere> spheres;
    Vec3 backgroundColor;

    Scene() : backgroundColor(0.2f, 0.2f, 0.2f) {
        // 添加球体
        spheres.push_back(Sphere(Vec3(0, 0, -5), 1.0, Material(DIFFUSE, Vec3(0.8f, 0.3f, 0.3f))));
        spheres.push_back(Sphere(Vec3(-2.5f, 0, -6), 1.0, Material(METAL, Vec3(0.8f, 0.8f, 0.8f), 0.1f)));
        spheres.push_back(Sphere(Vec3(2.5f, 0, -6), 1.0, Material(GLASS, Vec3(1.0f, 1.0f, 1.0f), 0.0f, 1.5f)));
        spheres.push_back(Sphere(Vec3(0, -500.5f, -5), 500.0, Material(DIFFUSE, Vec3(0.2f, 0.8f, 0.2f)))); // 地面
    }

    // 遍历场景求交
    bool intersect(const Ray& ray, float& t, int& id) const {
        float min_t = std::numeric_limits<float>::max();
        bool hit = false;

        for (int i = 0; i < spheres.size(); i++) {
            float d;
            if (spheres[i].intersect(ray, d) && d < min_t) {
                min_t = d;
                id = i;
                hit = true;
            }
        }

        t = min_t;
        return hit;
    }
};

// 反射计算
Vec3 reflect(const Vec3& v, const Vec3& n) {
    return v - n * (2.0f * v.dot(n));
}

// 折射计算
Vec3 refract(const Vec3& v, const Vec3& n, float ni_over_nt, bool& total_reflection) {
    Vec3 uv = v.normalize();
    float dt = uv.dot(n);
    float discriminant = 1.0f - ni_over_nt * ni_over_nt * (1 - dt * dt);
    if (discriminant > 0) {
        total_reflection = false;
        return (uv - n * dt) * ni_over_nt - n * std::sqrt(discriminant);
    }
    else {
        total_reflection = true;
        return Vec3(0, 0, 0);
    }
}

// 随机数生成
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<float> dis(0.0f, 1.0f);

float randomFloat() {
    return dis(gen);
}

// 计算光线颜色
Vec3 trace(const Ray& ray, const Scene& scene, int depth) {
    if (depth >= MAX_DEPTH) return Vec3(0, 0, 0); // 递归终止

    float t; // 交点距离
    int id;  // 最近物体ID

    if (scene.intersect(ray, t, id)) {
        const Sphere& obj = scene.spheres[id];
        Vec3 hitPoint = ray.pointAt(t);
        Vec3 normal = (hitPoint - obj.center).normalize();

        // 计算阴影偏移点
        Vec3 offsetPoint = hitPoint + normal * 0.001f;

        switch (obj.material.type) {
        case DIFFUSE: {
            // 简单环境光
            Vec3 color = obj.material.albedo * 0.2f;

            // 添加一个点光源
            Vec3 lightPos(-5, 5, -3);
            Vec3 toLight = (lightPos - offsetPoint).normalize();

            // 检查光线是否被阻挡
            Ray shadowRay(offsetPoint, toLight);
            float lightDist = (lightPos - offsetPoint).length();
            bool blocked = false;
            for (const Sphere& s : scene.spheres) {
                float tempT;
                if (s.intersect(shadowRay, tempT) && tempT < lightDist) {
                    blocked = true;
                    break;
                }
            }

            if (!blocked) {
                float diff = std::max(0.0f, normal.dot(toLight));
                color = color + obj.material.albedo * diff * 0.8f;
            }

            return color;
        }

        case METAL: {
            Ray reflected(offsetPoint, reflect(ray.direction, normal));
            // 添加随机偏移模拟粗糙度
            Vec3 jitter = Vec3(
                (randomFloat() - 0.5f) * 2.0f,
                (randomFloat() - 0.5f) * 2.0f,
                (randomFloat() - 0.5f) * 2.0f
            ).normalize() * obj.material.roughness;
            reflected.direction = (reflected.direction + jitter).normalize();
            return trace(reflected, scene, depth + 1).mult(obj.material.albedo);
        }

        case GLASS: {
            Vec3 outwardNormal;
            float ni_over_nt;
            Vec3 attenuation(1, 1, 1);
            Ray scattered(offsetPoint, Vec3(0, 0, 0));

            if (ray.direction.dot(normal) < 0) { // 外部光线
                outwardNormal = normal;
                ni_over_nt = 1.0f / obj.material.ior;
            }
            else { // 内部光线
                outwardNormal = normal * -1.0f;
                ni_over_nt = obj.material.ior;
            }

            bool total_reflection;
            Vec3 refracted = refract(ray.direction, outwardNormal, ni_over_nt, total_reflection);

            if (!total_reflection) {
                scattered = Ray(offsetPoint, refracted);
            }
            else {
                scattered = Ray(offsetPoint, reflect(ray.direction, normal));
            }

            return trace(scattered, scene, depth + 1).mult(attenuation);
        }
        }
    }

    return scene.backgroundColor; // 背景色
}

// 自定义的clamp函数
float clamp(float value, float min, float max) {
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

int main() {
    // 初始化SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "SDL初始化失败: " << SDL_GetError() << std::endl;
        return 1;
    }

    SDL_Window* window = SDL_CreateWindow("SDL3光线追踪演示",
        WIDTH, HEIGHT, SDL_WINDOW_RESIZABLE);

    if (!window) {
        std::cerr << "窗口创建失败: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return 1;
    }

    // 使用SDL3的标志创建渲染器
    SDL_Renderer* renderer = SDL_CreateRenderer(window,
        nullptr);

    if (!renderer) {
        std::cerr << "渲染器创建失败: " << SDL_GetError() << std::endl;
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    // 创建纹理
    SDL_Texture* texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888,
        SDL_TEXTUREACCESS_STREAMING, WIDTH, HEIGHT);
    if (!texture) {
        std::cerr << "纹理创建失败: " << SDL_GetError() << std::endl;
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    // 创建场景和相机
    Scene scene;
    Camera camera(Vec3(0, 0, 0), Vec3(0, 0, -1), 70.0f);

    // 像素缓冲区
    std::vector<Uint32> pixels(WIDTH * HEIGHT);

    // 渲染场景
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            // 将像素坐标转换为 -1到1 的范围
            float u = (x + 0.5f) / WIDTH * 2.0f - 1.0f;
            float v = (HEIGHT - y - 0.5f) / HEIGHT * 2.0f - 1.0f; // 翻转Y轴

            Ray ray = camera.getRay(u, v);
            Vec3 color = trace(ray, scene, 0);

            // 使用自定义clamp函数确保颜色值在有效范围内
            Uint32 r = static_cast<Uint32>(255.0f * clamp(color.x, 0.0f, 1.0f));
            Uint32 g = static_cast<Uint32>(255.0f * clamp(color.y, 0.0f, 1.0f));
            Uint32 b = static_cast<Uint32>(255.0f * clamp(color.z, 0.0f, 1.0f));

            pixels[y * WIDTH + x] = (0xFF << 24) | (r << 16) | (g << 8) | b;
        }
    }

    // 更新纹理
    SDL_UpdateTexture(texture, nullptr, pixels.data(), WIDTH * sizeof(Uint32));

    // 渲染循环
    bool running = true;
    while (running) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_EVENT_QUIT) {
                running = false;
            }
            if (event.type == SDL_EVENT_KEY_DOWN) {
                running = false;
            }
        }

        // 渲染到屏幕
        SDL_RenderClear(renderer);
        SDL_RenderTexture(renderer, texture, nullptr, nullptr);
        SDL_RenderPresent(renderer);

        SDL_Delay(16); // 约60FPS
    }

    // 清理资源
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}

//#include<SDL.h>
//
//int main(int argc, char* argv[])
//{
//	//初始化SDL
//	if (SDL_Init(SDL_INIT_VIDEO) < 0)
//	{
//		SDL_Log("can not init SDL:%s", SDL_GetError());
//		return -1;
//	}
//
//	return 0;
//}