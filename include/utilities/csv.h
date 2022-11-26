#ifndef csv_file_local_h
#define csv_file_local_h
#include <iostream>

using namespace std;
#include <fstream>
#include <stdio.h>

class csvfile
{
    std::ofstream fs_;
    const std::string separator_;
    const std::string filename_;
public:
    csvfile(const std::string filename, const std::string separator = ";")
        : fs_()
        , separator_(separator)
        , filename_(filename)
    {
        std::remove(filename_.c_str());
        fs_.exceptions(std::ios::failbit | std::ios::badbit);
        fs_.open(filename,  std::ios::trunc);
    }

    ~csvfile()
    {
        flush();
        fs_.close();
    }

    void flush()
    {
        fs_.flush();
    }

    void endrow()
    {
        fs_ << std::endl;
    }

    csvfile& operator << ( csvfile& (* val)(csvfile&))
    {
        return val(*this);
    }

    csvfile& operator << (const char * val)
    {
        fs_ << '"' << val << '"' << separator_;
        return *this;
    }

    csvfile& operator << (const std::string & val)
    {
        fs_ << '"' << val << '"' << separator_;
        return *this;
    }

    template<typename T>
    csvfile& operator << (const T& val)
    {
        fs_ << val << separator_;
        return *this;
    }
    
    void AddRow() 
    {
        fs_<<std::endl;
    }
    template<typename First, typename ... Strings>
    void AddRow(First arg, const Strings&... rest) 
    {
        fs_<<arg<<separator_;
        AddRow(rest...);
    }
};
#endif
