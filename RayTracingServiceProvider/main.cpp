#include <tuple>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <chrono>

struct vec3 {
    float x = 0, y = 0, z = 0;
    vec3() = default;
    vec3(float x, float y, float z) : x(x), y(y), z(z) {}

    float& operator[](int i) { return (&x)[i]; }
    const float& operator[](int i) const { return (&x)[i]; }

    vec3 operator*(float v) const { return { x * v, y * v, z * v }; }
    vec3 operator*(const vec3& v) const { return { x * v.x, y * v.y, z * v.z }; }
    vec3 operator+(const vec3& v) const { return { x + v.x, y + v.y, z + v.z }; }
    vec3 operator-(const vec3& v) const { return { x - v.x, y - v.y, z - v.z }; }
    vec3 operator-() const { return { -x, -y, -z }; }

    float dot(const vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    float norm() const { return std::sqrt(x * x + y * y + z * z); }
    vec3 normalized() const {
        float n = norm();
        return n > 0 ? (*this) * (1.f / n) : *this;
    }
};

vec3 operator*(float s, const vec3& v) {
    return v * s;
}

vec3 cross(const vec3& v1, const vec3& v2) {
    return {
        v1.y * v2.z - v1.z * v2.y,
        v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x
    };
}

struct Material {
    float refractive_index = 1;
    float albedo[4] = { 2, 0, 0, 0 };
    vec3 diffuse_color = { 0, 0, 0 };
    float specular_exponent = 0;
};

struct Cube {
    vec3 min_bound;  // 立方体最小顶点
    vec3 max_bound;  // 立方体最大顶点
    Material material;
};

const Material cube_material = { 1.0, {1.0f, 0.2f, 0.1f, 0.0f}, {0.2f, 0.7f, 0.2f}, 50.0f };

// 场景对象 (移除了所有球体，只保留一个立方体)
const std::vector<Cube> cubes = {
    {{0.5f, -1.5f, -15.0f}, {3.0f, 0.5f, -12.0f}, cube_material}  // 绿色立方体
};

const std::vector<vec3> lights = {
    {-20.0f, 20.0f, 20.0f},
    {30.0f, 50.0f, -25.0f},
    {30.0f, 20.0f, 30.0f}
};

vec3 reflect(const vec3& I, const vec3& N) {
    return I - 2.0f * I.dot(N) * N;
}

vec3 refract(const vec3& I, const vec3& N, float eta_t, float eta_i = 1.0f) {
    float cosi = -std::max(-1.0f, std::min(1.0f, I.dot(N)));
    if (cosi < 0) return refract(I, -N, eta_i, eta_t);

    float eta = eta_i / eta_t;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? vec3(1, 0, 0) : I * eta + N * (eta * cosi - std::sqrt(k));
}

// 立方体求交函数
std::pair<bool, float> ray_cube_intersect(const vec3& orig, const vec3& dir, const Cube& c) {
    // 计算射线进入和离开立方体的时间
    float t_min = -1e10f;
    float t_max = 1e10f;
    vec3 bounds[2] = { c.min_bound, c.max_bound };

    for (int i = 0; i < 3; i++) {
        float invD = 1.0f / dir[i];
        float t0 = (bounds[0][i] - orig[i]) * invD;
        float t1 = (bounds[1][i] - orig[i]) * invD;
        if (invD < 0.0f) std::swap(t0, t1);
        t_min = std::max(t_min, t0);  // 进入时间取最大值
        t_max = std::min(t_max, t1);  // 离开时间取最小值
        if (t_max < t_min) return { false, 0.0f };
    }

    // 找到相交点
    if (t_min < 0.0f) {
        if (t_max < 0.0f) return { false, 0.0f };
        return { true, t_max };
    }

    return { true, t_min };
}

struct IntersectionResult {
    bool hit;
    vec3 point;
    vec3 normal;
    Material material;
};

// 计算立方体相交点的法线
vec3 calculate_cube_normal(const vec3& point, const Cube& c) {
    const float epsilon = 0.0001f;

    // 检查点在立方体的哪个面上
    if (std::abs(point.x - c.min_bound.x) < epsilon) return { -1.0f, 0.0f, 0.0f };
    if (std::abs(point.x - c.max_bound.x) < epsilon) return { 1.0f, 0.0f, 0.0f };
    if (std::abs(point.y - c.min_bound.y) < epsilon) return { 0.0f, -1.0f, 0.0f };
    if (std::abs(point.y - c.max_bound.y) < epsilon) return { 0.0f, 1.0f, 0.0f };
    if (std::abs(point.z - c.min_bound.z) < epsilon) return { 0.0f, 0.0f, -1.0f };
    if (std::abs(point.z - c.max_bound.z) < epsilon) return { 0.0f, 0.0f, 1.0f };

    // 默认返回上表面法线
    return { 0.0f, 1.0f, 0.0f };
}

IntersectionResult scene_intersect(const vec3& orig, const vec3& dir) {
    vec3 pt, N;
    Material material;
    float nearest_dist = 1e10f;
    bool hit = false;

    // 地面检测
    if (std::abs(dir.y) > 1e-3f) {
        float d = -(orig.y + 4.0f) / dir.y;
        vec3 p = orig + dir * d;
        if (d > 1e-3f && d < nearest_dist && std::abs(p.x) < 10 && p.z < -10 && p.z > -30) {
            nearest_dist = d;
            pt = p;
            N = vec3(0, 1, 0);
            material.diffuse_color = (static_cast<int>(0.5f * pt.x + 1000) +
                static_cast<int>(0.5f * pt.z)) & 1 ?
                vec3(0.3f, 0.3f, 0.3f) : vec3(0.3f, 0.2f, 0.1f);
            hit = true;
        }
    }

    // 立方体检测 (唯一保留的物体)
    for (const auto& cube : cubes) {
        auto [intersects, dist] = ray_cube_intersect(orig, dir, cube);
        if (intersects && dist < nearest_dist) {
            nearest_dist = dist;
            pt = orig + dir * dist;
            N = calculate_cube_normal(pt, cube);
            material = cube.material;
            hit = true;
        }
    }

    return { hit, pt, N, material };
}

vec3 cast_ray(const vec3& orig, const vec3& dir, int depth = 0) {
    if (depth > 4) return vec3(0.2f, 0.7f, 0.8f); // 背景色

    auto intersection = scene_intersect(orig, dir);
    if (!intersection.hit) {
        return vec3(0.2f, 0.7f, 0.8f); // 背景色
    }

    vec3 point = intersection.point;
    vec3 N = intersection.normal;
    Material material = intersection.material;

    vec3 reflect_dir = reflect(dir, N).normalized();
    vec3 reflect_color = cast_ray(point + N * 1e-3f, reflect_dir, depth + 1);

    vec3 refract_dir = refract(dir, N, material.refractive_index).normalized();
    vec3 refract_color = cast_ray(point - N * 1e-3f, refract_dir, depth + 1);

    float diffuse = 0.0f;
    float specular = 0.0f;

    for (const auto& light : lights) {
        vec3 light_dir = (light - point).normalized();
        vec3 shadow_orig = point + N * 1e-3f;

        auto shadow_int = scene_intersect(shadow_orig, light_dir);
        if (shadow_int.hit &&
            (shadow_int.point - shadow_orig).norm() < (light - shadow_orig).norm()) {
            continue; // 在阴影中
        }

        diffuse += std::max(0.0f, light_dir.dot(N));

        vec3 reflect_dir = reflect(-light_dir, N);
        float cos_alpha = std::max(0.0f, -reflect_dir.dot(dir));
        specular += std::pow(cos_alpha, material.specular_exponent);
    }

    vec3 result =
        material.diffuse_color * diffuse * material.albedo[0] +
        vec3(1.0f, 1.0f, 1.0f) * specular * material.albedo[1] +
        reflect_color * material.albedo[2] +
        refract_color * material.albedo[3];

    // 限制颜色范围
    for (int i = 0; i < 3; i++) {
        if (result[i] < 0) result[i] = 0;
        if (result[i] > 1) result[i] = 1;
    }

    return result;
}

int main() {
    const int width = 1024;
    const int height = 768;
    const float fov = 1.05f; // 60度视野
    std::vector<vec3> framebuffer(width * height);

    auto start = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float dir_x = (x + 0.5f) - width / 2.0f;
            float dir_y = -(y + 0.5f) + height / 2.0f;
            float dir_z = -height / (2.0f * std::tan(fov / 2.0f));
            framebuffer[y * width + x] = cast_ray(vec3(0, 0, 0), vec3(dir_x, dir_y, dir_z).normalized());
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "渲染时间: " << duration.count() << " ms" << std::endl;

    std::ofstream ofs("outputcube.ppm", std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";

    for (auto& color : framebuffer) {
        // Gamma校正
        color.x = std::sqrt(color.x);
        color.y = std::sqrt(color.y);
        color.z = std::sqrt(color.z);

        ofs << static_cast<char>(255 * std::min(1.0f, color.x))
            << static_cast<char>(255 * std::min(1.0f, color.y))
            << static_cast<char>(255 * std::min(1.0f, color.z));
    }

    return 0;
}