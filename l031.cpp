#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cctype>
#include <random>
#include <string>
#include <iterator>
#include <list>
#include <algorithm>
#include <chrono>

using namespace std;

const int SIZE = 100000;
namespace ppm
{
    class Image
    {
    public:
        Image() {}
        Image(string name) : width(800), height(800), name(name)
        {
            imageData = new int[800][800 * 3];
            for (int i = 0; i < width; i++)
            {
                for (int j = 0; j < height * 3; j++)
                {
                    imageData[i][j] = 0;
                }
            }
        }
        void fill(int r, int g, int b)
        {
            for (int i = 0; i < width; i++)
            {
                for (int j = 0; j < height; j++)
                {
                    imageData[i][j * 3 + 0] = r;
                    imageData[i][j * 3 + 1] = g;
                    imageData[i][j * 3 + 2] = b;
                }
            }
        }
        void output()
        {
            // Create and open a text file
            ofstream out(this->name);

            // Write header to file
            out << "P3 " << this->width << " " << this->height << " 255\r\n";

            // Write data
            for (int i = 0; i < this->width; i++)
            {
                for (int j = 0; j < this->height * 3; j++)
                {
                    out << this->imageData[i][j] << " ";
                }
                out << "\n";
            }

            // Close the file
            out.close();
        }
        void write(int x, int y, int r, int g, int b)
        {
            if (x < 0 || y < 0 || x >= this->width || y >= this->height)
            {
                return;
            }

            this->imageData[y][x * 3 + 0] = r;
            this->imageData[y][x * 3 + 1] = g;
            this->imageData[y][x * 3 + 2] = b;
        }
        int get_width()
        {
            return this->width;
        }
        int get_height()
        {
            return this->height;
        }

    private:
        int width;
        int height;
        std::string name;
        int (*imageData)[800 * 3];
    };

    class Point
    {
    public:
        Point() {}

        Point(double x, double y) : x(x), y(y) {}

        double xpos()
        {
            return this->x;
        }

        double ypos()
        {
            return this->y;
        }

        void set_xpos(double x)
        {
            this->x = x;
        }

        void set_ypos(double y)
        {
            this->y = y;
        }
        double calc_dist(Point &other)
        {
            double dx, dy, dist;
            dx = this->x - other.xpos();
            dy = this->y - other.ypos();
            dist = sqrt(dx * dx + dy * dy);
            return dist;
        }

    private:
        double x, y;
    };

    class Circle
    {
    public:
        Circle(double x, double y, double r) : x(x), y(y), r(r)
        {
        }
        void draw(Image &image, int r, int g, int b)
        {
            int x, y, y2, y2_new, ty, xoff, yoff;
            x = 0;
            y = (int)(this->r * (double)image.get_width());
            y2 = y * y;
            y2_new = y2;
            ty = (2 * y) - 1;
            xoff = (this->x * (double)image.get_width());
            yoff = (this->y * (double)image.get_height());

            while (x <= y)
            {
                if ((y2 - y2_new) >= ty)
                {
                    y2 -= ty;
                    y -= 1;
                    ty -= 2;
                }
                image.write(xoff + x, yoff + y, r, g, b);
                image.write(xoff + x, yoff - y, r, g, b);
                image.write(xoff - x, yoff + y, r, g, b);
                image.write(xoff - x, yoff - y, r, g, b);
                image.write(xoff + y, yoff + x, r, g, b);
                image.write(xoff + y, yoff - x, r, g, b);
                image.write(xoff - y, yoff + x, r, g, b);
                image.write(xoff - y, yoff - x, r, g, b);
                y2_new -= (2 * x) - 3;
                x += 1;
            }
        }
        double radius()
        {
            return this->r;
        }
        double x_pos()
        {
            return this->x;
        }
        double y_pos()
        {
            return this->y;
        }

    private:
        double x, y, r;
    };
}

namespace gen
{
    using namespace ppm;

    bool get_answer()
    {
        while (true)
        {
            string s;
            cout << "Do you want to generate the points? (Yes/No)" << endl;
            cin >> s;
            auto first = toupper(s[0]);
            if (first == 'Y')
            {
                return true;
            }
            else if (first == 'N')
            {
                return false;
            }
            else
            {
                cout << "Invalid input, try again:" << endl;
            }
        }
    }

    void part0()
    {
        // Maybe Gen Points
        auto gen = get_answer();
        if (!gen)
        {
            return;
        }
        else
        {
            // Gen Points
            random_device os_seed;
            mt19937 gen(os_seed());
            uniform_real_distribution<> xy(0, 1);
            auto points = list<Point>(SIZE, Point());
            for (auto &point : points)
            {
                point.set_xpos(xy(gen));
                point.set_ypos(xy(gen));
            }
            // Output to a file
            ofstream output("points.txt");
            auto output_point = [&](Point &point)
            { output << point.xpos() << "  " << point.ypos() << endl; };
            output << fixed << setprecision(23);
            for (auto &point : points)
            {
                output_point(point);
            }
            output.close();
        }
    }
}

namespace compare
{
    using namespace ppm;

    list<Point> read_file()
    {
        list<Point> ret;
        ifstream file("points.txt");

        double x, y;
        while (file >> x >> y)
        {
            ret.push_front(Point(x, y));
        }

        return ret;
    }

    void find_min_dist(list<Point> &points, Point **p1, Point **p2)
    {
        auto min_dist = 2.0;

        for (auto it = points.begin(); it != points.end(); it++)
        {
            // Compare
            for (auto it2 = next(it); it2 != points.end(); it2++)
            {
                auto dist = it->calc_dist(*it2);
                if (dist < min_dist)
                {
                    min_dist = dist;
                    *p1 = &*it;
                    *p2 = &*it2;
                }
            }
        }
    }

    long long part1()
    {
        // Read points
        auto points = read_file();
        // Start timing
        auto start = chrono::high_resolution_clock::now();
        // Find closest
        Point *p1, *p2;
        find_min_dist(points, &p1, &p2);
        // End timing
        auto elapsed = chrono::high_resolution_clock::now() - start;
        // Get length
        long long microseconds = chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        // Draw them and output points
        ofstream results("results.txt");
        results << fixed << setprecision(10) << "Found Points: (" << p1->xpos() << ", " << p1->ypos() << "), (";
        results << p2->xpos() << ", " << p2->ypos() << ")" << endl;
        results << "They have a distance of " << p1->calc_dist(*p2) << endl;
        results << "It took " << microseconds << " microseconds to calculate" << endl;
        Image output("points.ppm");
        output.fill(255, 255, 255);
        for (auto &point : points)
        {
            if (&point != p1 && &point != p2)
            {
                Circle(point.xpos(), point.ypos(), 2.0 / 800.0).draw(output, 0, 0, 0);
                Circle(point.xpos(), point.ypos(), 3.0 / 800.0).draw(output, 0, 0, 0);
            }
            else
            {
                Circle(point.xpos(), point.ypos(), 2.0 / 800.0).draw(output, 255, 0, 0);
                Circle(point.xpos(), point.ypos(), 3.0 / 800.0).draw(output, 255, 0, 0);
            }
        }
        output.output();
        return microseconds;
    }

    void split(Point *points, int len, Point **p1, Point **p2, double *min_dist)
    {
        // Ending Condition
        if (len <= 1)
        {
            // There's only 1 point, so this shouldn't ever happen
            int a = 5;
        }
        else if (len <= 2)
        {
            // There's 2 points
            *min_dist = points[0].calc_dist(points[1]);
            *p1 = &points[0];
            *p2 = &points[1];
        }
        else if (len <= 3)
        {
            // There's 3 points
            double d1, d2, d3;
            d1 = points[0].calc_dist(points[1]);
            d2 = points[1].calc_dist(points[2]);
            d3 = points[2].calc_dist(points[0]);
            // Find min
            if (d1 < d2 && d1 < d3)
            {
                *min_dist = d1;
                *p1 = &points[0];
                *p2 = &points[1];
            }
            else if (d2 < d3)
            {
                *min_dist = d2;
                *p1 = &points[1];
                *p2 = &points[2];
            }
            else
            {
                *min_dist = d3;
                *p1 = &points[2];
                *p2 = &points[0];
            }
        }
        else
        {
            // Split Left
            int len2 = len / 2;
            split(&points[0], len2, p1, p2, min_dist);
            // Split right
            Point *p3, *p4;
            double min_dist2;
            split(&points[len2], len - len2, &p3, &p4, &min_dist2);
            // Get min
            if (*min_dist > min_dist2)
            {
                *min_dist = min_dist2;
                *p1 = p3;
                *p2 = p4;
            }
            // Combine middle
            int i = len2;
            while (i >= 0 && points[i].xpos() > points[len2].xpos() - *min_dist)
            {
                int j = len2;
                if (i == j)
                {
                    j++;
                }
                while (j < len && points[j].xpos() < points[len2].xpos() + *min_dist)
                {
                    double dist = points[i].calc_dist(points[j]);
                    if (dist < *min_dist)
                    {
                        *min_dist = dist;
                        *p1 = &points[i];
                        *p2 = &points[j];
                    }
                    j++;
                }
                i--;
            }
            // Done
        }
    }

    long long part2()
    {
        // Read points
        auto points_list = read_file();
        // Turn into a vector
        auto points = vector<Point>(points_list.begin(), points_list.end());
        // Start timing
        auto start = chrono::high_resolution_clock::now();
        // Sort by x coord
        auto xmax = [](Point &p1, Point &p2)
        { return p1.xpos() < p2.xpos(); };
        sort(points.begin(), points.end(), xmax);
        // Recurse
        Point *p1, *p2;
        double min_dist;
        split(&points[0], points.size(), &p1, &p2, &min_dist);
        // End timing
        auto elapsed = chrono::high_resolution_clock::now() - start;
        // Get length
        long long microseconds = chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        // Draw them and output points
        ofstream results("results2.txt");
        results << fixed << setprecision(10) << "Found Points: (" << p1->xpos() << ", " << p1->ypos() << "), (";
        results << p2->xpos() << ", " << p2->ypos() << ")" << endl;
        results << "They have a distance of " << p1->calc_dist(*p2) << endl;
        results << "It took " << microseconds << " microseconds to calculate" << endl;
        Image output("points2.ppm");
        output.fill(255, 255, 255);
        for (auto &point : points)
        {
            if (&point != p1 && &point != p2)
            {
                Circle(point.xpos(), point.ypos(), 2.0 / 800.0).draw(output, 0, 0, 0);
                Circle(point.xpos(), point.ypos(), 3.0 / 800.0).draw(output, 0, 0, 0);
            }
            else
            {
                Circle(point.xpos(), point.ypos(), 2.0 / 800.0).draw(output, 255, 0, 0);
                Circle(point.xpos(), point.ypos(), 3.0 / 800.0).draw(output, 255, 0, 0);
            }
        }
        output.output();
        return microseconds;
    }
}

int main()
{
    gen::part0();

    auto len1 = compare::part1();
    cout << "Part 1 took " << len1 << " microseconds" << endl
         << endl;

    auto len2 = compare::part2();
    cout << "Part 2 took " << len2 << " microseconds" << endl;

    return 0;
}