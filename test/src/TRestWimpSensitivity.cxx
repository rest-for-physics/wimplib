
#include <TRestWimpSensitivity.h>
#include <gtest/gtest.h>

#include <filesystem>

namespace fs = std::filesystem;

using namespace std;

const auto filesPath = fs::path(__FILE__).parent_path().parent_path() / "files";
const auto wimpRml = filesPath / "wimp.rml";

TEST(TRestWimpSensitivity, TestFiles) {
    cout << "Test files path: " << filesPath << endl;

    // Check dir exists and is a directory
    EXPECT_TRUE(fs::is_directory(filesPath));
    // Check it's not empty
    EXPECT_TRUE(!fs::is_empty(filesPath));

    // All used files in this tests
    EXPECT_TRUE(fs::exists(wimpRml));
}

TEST(TRestWimpSensitivity, FromRml) {
    cout << "Path: " << wimpRml << endl;

    TRestWimpSensitivity WS(wimpRml.c_str());

    WS.PrintMetadata();

    EXPECT_TRUE(abs(WS.GetSensitivity(1) - 2.9732007e-39) < 1E-45 );
}
