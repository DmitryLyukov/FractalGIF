#include <string>
#include <cstdint>
#include <vector>
#include <cmath>
#include <complex>
#include <iostream>

#include "CImg.h"

#include <thread>
#include <atomic>

int iter_for_point(const std::complex<double> &c,
        const float power, const int MAX_ITER = 500) {
    
    int i = 0;
    std::complex<double> z_i = 0;
    double sqr_abs = pow(real(z_i), 2) + pow(imag(z_i), 2);

    while ((i < MAX_ITER) && (sqr_abs < 4)) {
        std::complex<double> sqr = pow(z_i, power);
        z_i = sqr + c;
        ++i;
        sqr_abs = pow(real(z_i), 2) + pow(imag(z_i), 2);
    }
    
    return i;
}


double norm_iter_for_point(const double Re, const double Im,
        const float power, const int MAX_ITER) {
    
    std::complex<double> z(Re, Im);
    return static_cast<double>(iter_for_point(z, power, MAX_ITER)) / MAX_ITER;
}


void painting(const double x1, const double x2,
        const double y1, const double y2,
        const int32_t width, const int32_t height,
        const char* file_name, const float power) {
    
    cimg_library::CImg<unsigned char> img(height, width, 1, 3);
    
    const double dy = (y2 - y1) / height;
    const double dx = (x2 - x1) / width;
    
    double y = y2;
    
    for (int32_t i = 0; i < height / 2; ++i) {
        double x = x1;
        
        for (int32_t j = 0; j < width; ++j) {
            const double norm_pnt = 1 - norm_iter_for_point(x, y, power);
            x += dx;
            
            const unsigned char red   = static_cast<unsigned char>(pow(sin(4.6 * norm_pnt), 2) * 255);
            const unsigned char green = static_cast<unsigned char>(pow(sin(4.6 * norm_pnt * norm_pnt * norm_pnt), 2) * 255);
            const unsigned char blue  = static_cast<unsigned char>(pow(sin(4.6 * norm_pnt), 2) * 255);
            
            img(i, j, 0) = red;
            img(i, j, 1) = green;
            img(i, j, 2) = blue;
            
            img(height - i - 1, j, 0) = red;
            img(height - i - 1, j, 1) = green;
            img(height - i - 1, j, 2) = blue;
        }
        
        y -= dy;
    }
    img.save_png(file_name);
}


void painting_thread(std::atomic<int64_t> &atom_cnt) {
    while (atom_cnt < 1605) {
        int64_t j = atom_cnt.fetch_add(1);
        if (j >= 1605) {
            break;
        }
        float i = 1. + 0.0025 * j;
        
        std::string str = std::to_string(j);
        while (str.length() < 5) {
            str = '0' + str;
        }
        str = str + '_' + std::to_string(j) + ".png";
        painting(-2, 1.2, -2.84, 2.84, 2160, 3840, str.c_str(), i);
        std::cout << j << std::endl;
    }
}

int main(int argc, char* argv[]) {
        
    size_t num_of_threads;
    if (argc > 1) {
        num_of_threads = atol(argv[1]);
    } else {
        num_of_threads = std::thread::hardware_concurrency();
    }
    if (num_of_threads == 0) {
        num_of_threads = 1;
    }

    std::atomic<int64_t> atom_cnt(0);
    
    std::vector<std::thread> threads;
    for (size_t i = 0; i < num_of_threads; ++i) {
        threads.push_back(std::thread (painting_thread, std::ref(atom_cnt)));
    }
    
    for (auto& thr: threads) {
        thr.join();
    }
    
    return 0;
}
