//
// Test Suite for geos::index::chain::MonotoneChain class.

#include <tut/tut.hpp>
#include <utility.h>
// geos
#include <geos/index/chain/MonotoneChain.h>


namespace tut {
// dummy data, not used
struct test_monotonechain_data {
    geos::io::WKTReader wktreader;

    std::shared_ptr<LineString> upwardWaveGeom;
    std::shared_ptr<LineString> downwardWaveGeom;
    

    test_monotonechain_data() {

        // ConvexPolygon2 wave{{0.00, -10.00}, {1.00, -9.88}, {2.00, -9.51}, {3.00, -8.91}, {4.00, -8.09}, {5.00, -7.07}, {6.00, -5.88}, {7.00, -4.54}, {8.00, -3.09}, {9.00, -1.56}, {10.00, 0.00}, {11.00, 1.56}, {12.00, 3.09}, {13.00, 4.54}, {14.00, 5.88}, {15.00, 7.07}, {16.00, 8.09}, {17.00, 8.91}, {18.00, 9.51}, {19.00, 9.88}};
        upwardWaveGeom = wktreader.read<LineString>(
            "LINESTRING(0.00 -10.00, 1.00 -9.88, 2.00 -9.51, 3.00 -8.91, 4.00 -8.09, 5.00 -7.07, 6.00 -5.88, "
            "7.00 -4.54, 8.00 -3.09, 9.00 -1.56, 10.00 -0.00, 11.00 1.56, 12.00 3.09, 13.00 4.54, 14.00 5.88, "
            "15.00 7.07, 16.00 8.09, 17.00 8.91, 18.00 9.51, 19.00 9.88)");

        // ConvexPolygon2 wave{{0.00, 10.00}, {1.00, 9.88}, {2.00, 9.51}, {3.00, 8.91}, {4.00, 8.09}, {5.00, 7.07}, {6.00, 5.88}, {7.00, 4.54}, {8.00, 3.09}, {9.00, 1.56}, {10.00, 0.00}, {11.00, -1.56}, {12.00, -3.09}, {13.00, -4.54}, {14.00, -5.88}, {15.00, -7.07}, {16.00, -8.09}, {17.00, -8.91}, {18.00, -9.51}, {19.00, -9.88}};
        downwardWaveGeom = wktreader.read<LineString>(
            "LINESTRING(0.00 10.00, 1.00 9.88, 2.00 9.51, 3.00 8.91, 4.00 8.09, 5.00 7.07, 6.00 5.88, 7.00 4.54, "
            "8.00 3.09, 9.00 1.56, 10.00 0.00, 11.00 -1.56, 12.00 -3.09, 13.00 -4.54, 14.00 -5.88, 15.00 -7.07, "
            "16.00 -8.09, 17.00 -8.91, 18.00 -9.51, 19.00 -9.88)");
    }
};

typedef test_group<test_monotonechain_data> group;
typedef group::object object;

group test_monotonechain_group("geos::index::chain::MonotoneChain");

void testWithLineStrings(const LineString& ls0, const LineString& ls1, bool expectedResult) {

    geos::index::chain::MonotoneChain ms0(*ls0.getCoordinatesRO(), 0, ls0.getNumPoints() - 1, nullptr);
    geos::index::chain::MonotoneChain ms1(*ls1.getCoordinatesRO(), 0, ls1.getNumPoints() - 1, nullptr);

    std::shared_ptr<LineString> ls0Reversed = ls0.reverse();
    std::shared_ptr<LineString> ls1Reversed = ls1.reverse();
    geos::index::chain::MonotoneChain ms0Reversed(*ls0Reversed->getCoordinatesRO(), 0, ls0Reversed->getNumPoints() - 1, nullptr);
    geos::index::chain::MonotoneChain ms1Reversed(*ls1Reversed->getCoordinatesRO(), 0, ls1Reversed->getNumPoints() - 1, nullptr);

    ensure_equals(ms0.intersects(ms1), expectedResult);
    ensure_equals(ms0.intersects(ms1Reversed), expectedResult);
    ensure_equals(ms0Reversed.intersects(ms1), expectedResult);
    ensure_equals(ms0Reversed.intersects(ms1Reversed), expectedResult);

    ensure_equals(ms1.intersects(ms0), expectedResult);
    ensure_equals(ms1.intersects(ms0Reversed), expectedResult);
    ensure_equals(ms1Reversed.intersects(ms0), expectedResult);
    ensure_equals(ms1Reversed.intersects(ms0Reversed), expectedResult);
}

//
// Test Cases

template<>
template<>
void object::test<1>
()
{
    std::shared_ptr<LineString> geom = wktreader.read<LineString>("LINESTRING(14.66 3.45, 17.61 12.63)");
    testWithLineStrings(*upwardWaveGeom, *geom, true);
}

template<>
template<>
void object::test<2>
()
{
    std::shared_ptr<LineString> geom = wktreader.read<LineString>("LINESTRING(2.15 -3.42, 9.51 11.36)");
    testWithLineStrings(*downwardWaveGeom, *geom, true);

}

template<>
template<>
void object::test<3>
()
{
    // ConvexPolygon2 geom{{7.62, -4.29}, {10.74, 0.38}};
    std::shared_ptr<LineString> geom = wktreader.read<LineString>("LINESTRING(7.62 -4.29, 10.74 0.38)");
    testWithLineStrings(*upwardWaveGeom, *geom, false);
}

template<>
template<>
void object::test<4>
()
{
    std::shared_ptr<LineString> geom = wktreader.read<LineString>("LINESTRING(6.26 6.98, 15.12 -6.12)");
    testWithLineStrings(*downwardWaveGeom, *geom, false);
}

template <>
template <>
void object::test<5>
()
{
    std::shared_ptr<LineString> geom = wktreader.read<LineString>("LINESTRING(4.23 -8.00, 4.74 -7.07)");
    testWithLineStrings(*upwardWaveGeom, *geom, true);
}

template<>
template<>
void object::test<6>
()
{
    std::shared_ptr<LineString> geom = wktreader.read<LineString>("LINESTRING(7.27 4.50, 7.67 3.29)");
    testWithLineStrings(*downwardWaveGeom, *geom, true);
}

template<>
template<>
void object::test<7>
()
{
    // ConvexPolygon2 polygon{{2.16, -11.13}, {5.51, -8.77}, {7.14, -5.64}, {8.51, -3.80}, {10.87, 0.19}, {12.37, 4.22}, {15.80, 9.06}};
    std::shared_ptr<LineString> geom = wktreader.read<LineString>(
        "LINESTRING(2.16 -11.13, 5.51 -8.77, 7.14 -5.64, 8.51 -3.80, 10.87 0.19, 12.37 4.22, 15.80 9.06)");
    testWithLineStrings(*upwardWaveGeom, *geom, true);
}

template<>
template<>
void object::test<8>
()
{
    // ConvexPolygon2 polygon{{2.45, 10.69}, {5.38, 8.55}, {6.92, 6.64}, {8.60, 3.80}, {10.41, 1.89}, {12.23, -4.40}, {14.61, -4.87}, {17.12, -6.55}};

    std::shared_ptr<LineString> geom = wktreader.read<LineString>(
        "LINESTRING(2.45 10.69, 5.38 8.55, 6.92 6.64, 8.60 3.80, 10.41 1.89, 12.23 -4.40, 14.61 -4.87, 17.12 -6.55)");
    testWithLineStrings(*downwardWaveGeom, *geom, true);
}

} // namespace tut
