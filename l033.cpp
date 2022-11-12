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
#include <unordered_map>
#include <utility>

using namespace std;

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
        string name;
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
        double calc_square_dist(Point &other)
        {
            double dx, dy;
            dx = this->x - other.xpos();
            dy = this->y - other.ypos();
            return dx * dx + dy * dy;
        }

        bool operator==(const Point &other)
        {
            return this->x == other.x && this->y == other.y;
        }

        bool operator!=(const Point &other)
        {
            return !(*this == other);
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

    bool get_answer(string question)
    {
        while (true)
        {
            string s;
            cout << question << endl;
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

    int get_size()
    {
        int size;
        cout << "How many? ";
        cin >> size;
        cout << endl;
        return size;
    }

    void gen_points(int size)
    {
        random_device os_seed;
        mt19937 gen(os_seed());
        uniform_real_distribution<> xy(0, 1);
        auto points = list<Point>(size, Point());
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

    bool part0()
    {
        // Ask if you want to output a PPM
        auto ppm = get_answer("Do you want to generate the PPM? (Yes/No)");
        // Maybe Gen Points
        auto gen = get_answer("Do you want to generate the points? (Yes/No)");
        if (!gen)
        {
            return ppm;
        }
        else
        {
            // Gen Points
            auto size = get_size();
            gen_points(size);
        }
        return ppm;
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

    template <class T>
    void find_min_dist(T &points, Point **p1, Point **p2)
    {
        auto min_dist = 2.0;

        for (auto it = points.begin(); it != points.end(); it++)
        {
            // Compare
            for (auto it2 = next(it); it2 != points.end(); it2++)
            {
                auto dist = it->calc_square_dist(*it2);
                if (dist < min_dist)
                {
                    min_dist = dist;
                    *p1 = &*it;
                    *p2 = &*it2;
                }
            }
        }
    }

    template <class T>
    long long time_part1(T &points, Point *point1, Point *point2)
    {
        // Start timing
        auto start = chrono::high_resolution_clock::now();
        // Find closest
        Point *p1 = nullptr, *p2 = nullptr;
        find_min_dist(points, &p1, &p2);
        *point1 = *p1;
        *point2 = *p2;
        // End timing
        auto elapsed = chrono::high_resolution_clock::now() - start;
        // Get length
        long long microseconds = chrono::duration_cast<chrono::microseconds>(elapsed).count();
        return microseconds;
    }

    long long part1(Point *point1, Point *point2, bool gen_ppm)
    {
        // Read points
        auto points = read_file();
        // Get time
        auto microseconds = time_part1(points, point1, point2);
        // Draw them
        if (gen_ppm)
        {
            Image output("points.ppm");
            output.fill(255, 255, 255);
            for (auto &point : points)
            {
                if (&point != point1 && &point != point2)
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
        }
        return microseconds;
    }

    void split(Point *points, int len, Point *p1, Point *p2, double *min_dist)
    {
        // Ending Condition
        if (len <= 1)
        {
            // There's only 1 point, so this shouldn't ever happen
        }
        else if (len <= 2)
        {
            // There's 2 points
            *min_dist = points[0].calc_square_dist(points[1]);
            *p1 = points[0];
            *p2 = points[1];
        }
        else if (len <= 3)
        {
            // There's 3 points
            double d1, d2, d3;
            d1 = points[0].calc_square_dist(points[1]);
            d2 = points[1].calc_square_dist(points[2]);
            d3 = points[2].calc_square_dist(points[0]);
            // Find min
            if (d1 < d2 && d1 < d3)
            {
                *min_dist = d1;
                *p1 = points[0];
                *p2 = points[1];
            }
            else if (d2 < d3)
            {
                *min_dist = d2;
                *p1 = points[1];
                *p2 = points[2];
            }
            else
            {
                *min_dist = d3;
                *p1 = points[2];
                *p2 = points[0];
            }
        }
        else
        {
            // Split Left
            int len2 = len / 2;
            split(&points[0], len2, p1, p2, min_dist);
            // Split right
            Point p3, p4;
            double min_dist2 = 0;
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
            while (i >= 0 && points[i].xpos() > points[len2].xpos() - sqrt(*min_dist))
            {
                int j = len2;
                if (i == j)
                {
                    j++;
                }
                while (j < len && points[j].xpos() < points[len2].xpos() + sqrt(*min_dist))
                {
                    double dist = points[i].calc_square_dist(points[j]);
                    if (dist < *min_dist)
                    {
                        *min_dist = dist;
                        *p1 = points[i];
                        *p2 = points[j];
                    }
                    j++;
                }
                i--;
            }
            // Done
        }
    }

    long long time_part2(vector<Point> &points, Point *point1, Point *point2)
    {
        // Start timing
        auto start = chrono::high_resolution_clock::now();
        // Sort by x coord
        auto xmax = [](Point &p1, Point &p2)
        { return p1.xpos() < p2.xpos(); };
        sort(points.begin(), points.end(), xmax);
        // Recurse
        double min_dist;
        split(&points[0], points.size(), point1, point2, &min_dist);
        // End timing
        auto elapsed = chrono::high_resolution_clock::now() - start;
        // Get length
        long long microseconds = chrono::duration_cast<chrono::microseconds>(elapsed).count();
        return microseconds;
    }

    long long part2(Point *point1, Point *point2, bool output_ppm)
    {
        // Read points
        auto points_list = read_file();
        // Turn into a vector
        auto points = vector<Point>(points_list.begin(), points_list.end());
        // Get time
        auto microseconds = time_part2(points, point1, point2);
        // Draw them
        if (output_ppm)
        {
            Image output("points2.ppm");
            output.fill(255, 255, 255);
            for (auto &point : points)
            {
                if (point != *point1 && point != *point2)
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
        }
        return microseconds;
    }

    void split3(Point *points, int len, Point *p1, Point *p2, double *min_dist)
    {
        // Ending Condition
        if (len <= 1)
        {
            // There's only 1 point, so this shouldn't ever happen
        }
        else if (len <= 2)
        {
            // There's 2 points
            *min_dist = points[0].calc_square_dist(points[1]);
            *p1 = points[0];
            *p2 = points[1];
        }
        else if (len <= 3)
        {
            // There's 3 points
            double d1, d2, d3;
            d1 = points[0].calc_square_dist(points[1]);
            d2 = points[1].calc_square_dist(points[2]);
            d3 = points[2].calc_square_dist(points[0]);
            // Find min
            if (d1 < d2 && d1 < d3)
            {
                *min_dist = d1;
                *p1 = points[0];
                *p2 = points[1];
            }
            else if (d2 < d3)
            {
                *min_dist = d2;
                *p1 = points[1];
                *p2 = points[2];
            }
            else
            {
                *min_dist = d3;
                *p1 = points[2];
                *p2 = points[0];
            }
        }
        else
        {
            // Split Left
            int len2 = len / 2;
            split3(&points[0], len2, p1, p2, min_dist);
            // Split right
            Point p3, p4;
            double min_dist2 = 0;
            split3(&points[len2], len - len2, &p3, &p4, &min_dist2);
            // Get min
            if (*min_dist > min_dist2)
            {
                *min_dist = min_dist2;
                *p1 = p3;
                *p2 = p4;
            }
            // Combine middle
            auto max = len2, min = len2 + 1;
            for (; min > 0 && points[len2].xpos() - points[min].xpos() < sqrt(*min_dist); min--)
            {
                for (; max < len && points[max].xpos() - points[len2].xpos() < sqrt(*min_dist); max++)
                {
                }
            }
            vector<Point> strip(&points[min], &points[max]);
            auto ymax = [](Point &p1, Point &p2)
            { return p1.ypos() < p2.ypos(); };
            sort(strip.begin(), strip.end(), ymax);
            for (unsigned int i = 0; i < strip.size(); i++)
            {
                for (unsigned int j = i + 1; j < strip.size() && j < i + 8; j++)
                {
                    auto dist = strip[i].calc_square_dist(strip[j]);
                    if (dist < *min_dist)
                    {
                        *min_dist = dist;
                        *p1 = strip[i];
                        *p2 = strip[j];
                    }
                }
            }
            // Done
        }
    }

    long long time_part3(vector<Point> &points, Point *point1, Point *point2)
    {
        // Start timing
        auto start = chrono::high_resolution_clock::now();
        // Sort by x coord
        auto xmax = [](Point &p1, Point &p2)
        { return p1.xpos() < p2.xpos(); };
        sort(points.begin(), points.end(), xmax);
        // Recurse
        double min_dist;
        split3(&points[0], points.size(), point1, point2, &min_dist);
        // End timing
        auto elapsed = chrono::high_resolution_clock::now() - start;
        // Get length
        long long microseconds = chrono::duration_cast<chrono::microseconds>(elapsed).count();
        return microseconds;
    }

    long long part3(Point *point1, Point *point2, bool output_ppm)
    {
        // Read points
        auto points_list = read_file();
        // Turn into a vector
        auto points = vector<Point>(points_list.begin(), points_list.end());
        // Get time
        auto microseconds = time_part3(points, point1, point2);
        // Draw them
        if (output_ppm)
        {
            Image output("points3.ppm");
            output.fill(255, 255, 255);
            for (auto &point : points)
            {
                if (point != *point1 && point != *point2)
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
        }
        return microseconds;
    }

    long long time_part4(vector<Point> &points, Point *p1, Point *p2)
    {
        // Start timing
        auto start = chrono::high_resolution_clock::now();
        // Make map with enough space for all points
        class Hash
        {
        public:
            size_t operator()(const pair<long long, long long> &key) const
            {
                auto first = hash<long long>{}(key.first);
                auto second = hash<long long>{}(key.second);
                if (first != second)
                {
                    return first ^ second;
                }
                return first;
            }
        };
        unordered_map<pair<long long, long long>, Point, Hash> map(points.size());
        // Set initial variables
        double min_dist = 1.0;
        auto insert_point = [&](Point point)
        { map[make_pair(point.xpos() * 2.0 / min_dist, point.ypos() * 2.0 / min_dist)] = point; };
        // Start iteration
        for (Point &point : points)
        {
            // Check surrounding squares
            long long i = point.xpos() * 2.0 / min_dist;
            long long j = point.ypos() * 2.0 / min_dist;
            double new_min_dist = 1.0;
            for (long long x = -2; x <= 2; x++)
            {
                for (long long y = -2; y <= 2; y++)
                {
                    if (map.find(make_pair(i + x, j + y)) != map.end())
                    {
                        auto p = map[make_pair(i + x, j + y)];
                        new_min_dist = point.calc_dist(p);
                        *p1 = point;
                        *p2 = p;
                    }
                }
            }
            // Maybe rehash stuff
            if (new_min_dist < min_dist)
            {
                min_dist = new_min_dist;
                map.clear();
                for (auto it = points.begin(); *it != point; it++)
                {
                    insert_point(*it);
                }
            }
            // Insert new point
            insert_point(point);
        }
        // End timing
        auto elapsed = chrono::high_resolution_clock::now() - start;
        // Get length
        long long microseconds = chrono::duration_cast<chrono::microseconds>(elapsed).count();
        return microseconds;
    }

    long long part4(Point *point1, Point *point2, bool output_ppm)
    {
        // Read points
        auto points_list = read_file();
        // Turn into a vector
        auto points = vector<Point>(points_list.begin(), points_list.end());
        // Get time
        auto microseconds = time_part4(points, point1, point2);
        // Draw them
        if (output_ppm)
        {
            Image output("points4.ppm");
            output.fill(255, 255, 255);
            for (auto &point : points)
            {
                if (point != *point1 && point != *point2)
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
        }
        return microseconds;
    }

    void gen_csv(bool part1, bool part2, bool part3, bool part4)
    {
        ofstream csv("results.csv");
        csv << "# Of Points";
        if (part1)
        {
            csv << ";part1";
        }
        if (part2)
        {
            csv << ";part2";
        }
        if (part3)
        {
            csv << ";part3";
        }
        if (part4)
        {
            csv << ";part4";
        }
        csv << endl
            << "0;0;0" << endl;
        random_device os_seed;
        mt19937 gen(os_seed());
        uniform_real_distribution<> xy(0, 1);
        auto num = 50000000;
        auto points_all = vector<Point>(num, Point());
        for (auto &point : points_all)
        {
            point.set_xpos(xy(gen));
            point.set_ypos(xy(gen));
        }
        auto incr = 1000000;
        for (int count = incr; count <= num; count += incr)
        {
            cout << count << " Points... " << endl;
            csv << count;
            /* Gen */
            Point p1, p2;
            auto points = vector<Point>(points_all.begin(), points_all.begin() + count);
            /* Part 1 */
            if (part1)
            {
                auto r1 = compare::time_part1(points, &p1, &p2);
                csv << ";" << r1;
            }
            /* Part 2 */
            if (part2)
            {
                auto r2 = compare::time_part2(points, &p1, &p2);
                csv << ";" << r2;
            }
            /* Part 3 */
            if (part3)
            {
                auto r3 = compare::time_part3(points, &p1, &p2);
                csv << ";" << r3;
            }
            /* Part 4 */
            if (part4)
            {
                auto r4 = compare::time_part4(points, &p1, &p2);
                csv << ";" << r4;
            }
            /* Finish Line */
            csv << endl;
        }
        csv.close();
    }
}

int main()
{
    compare::gen_csv(false, false, true, true);
    // bool output_ppm = gen::part0();

    // ofstream results("results.txt");
    // ppm::Point p1, p2;
    // char buf[230];

    // auto len1 = compare::part1(&p1, &p2);
    // sprintf(buf,
    //         "Part 1 Found Points: (%.20f, %.20f), (%.20f, %.20f)\n"
    //         "They have a distance of: %f\n"
    //         "It took %lld microseconds to calculate.\n",
    //         p1.xpos(), p1.ypos(), p2.xpos(), p2.ypos(), p1.calc_dist(p2), len1);
    // results << buf << endl;
    // cout << buf << endl;

    // p1 = ppm::Point(0.0, 0.0);
    // p2 = ppm::Point(0.0, 0.0);
    // auto len2 = compare::part2(&p1, &p2, output_ppm);
    // sprintf(buf,
    //         "Part 2 Found Points: (%.20f, %.20f), (%.20f, %.20f)\n"
    //         "They have a distance of: %e\n"
    //         "It took %lld microseconds to calculate.\n",
    //         p1.xpos(), p1.ypos(), p2.xpos(), p2.ypos(), p1.calc_dist(p2), len2);
    // results << buf << endl;
    // cout << buf << endl;

    // p1 = ppm::Point(0.0, 0.0);
    // p2 = ppm::Point(0.0, 0.0);
    // auto len3 = compare::part3(&p1, &p2, output_ppm);
    // sprintf(buf,
    //         "Part 3 Found Points: (%.20f, %.20f), (%.20f, %.20f)\n"
    //         "They have a distance of: %e\n"
    //         "It took %lld microseconds to calculate.\n",
    //         p1.xpos(), p1.ypos(), p2.xpos(), p2.ypos(), p1.calc_dist(p2), len3);
    // results << buf;
    // cout << buf;

    // p1 = ppm::Point(0.0, 0.0);
    // p2 = ppm::Point(0.0, 0.0);
    // auto len4 = compare::part4(&p1, &p2, output_ppm);
    // sprintf(buf,
    //         "Part 4 Found Points: (%.20f, %.20f), (%.20f, %.20f)\n"
    //         "They have a distance of: %e\n"
    //         "It took %lld microseconds to calculate.\n",
    //         p1.xpos(), p1.ypos(), p2.xpos(), p2.ypos(), p1.calc_dist(p2), len4);
    // results << buf;
    // cout << buf;

    return 0;
}