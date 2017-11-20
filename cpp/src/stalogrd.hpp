#ifndef __NGPT_IGS_LOG_FILE__
#define __NGPT_IGS_LOG_FILE__

#include <string>
#include <fstream>
#include "ggdatetime/datetime_write.hpp"
#include "ggdatetime/datetime_read.hpp"

namespace ngpt
{

namespace igs_log_details
{

int
str_is_empty(const char* str)
{
    const char *c = str;
    while (c) {
        if (*c != ' ') return false;
        ++c;
    }
    return true;
}

/// Resolve a datetime from a string of type: 'CCYY-MM-DDThh:mmZ'
/// or 'CCYY-MM-DD'
template<typename T>
    ngpt::datetime<T> strptime_log(const char* str, char** stop=nullptr)
{
    char *end;
    const char* start = str;
    int ints[5];
    // std::cout<<"\n\tthis is what i have: ["<<str<<"]";

    for (int i = 0; i < 5; ++i) {
        ints[i] = static_cast<int>( std::abs(std::strtol(start, &end, 10)) );
        // std::cout<<"\n\tend now set to: "<<*end;
        if (errno == ERANGE || start == end) {
            errno = 0;
            throw std::invalid_argument
                ("Invalid date format: \""+std::string(str)+"\" (argument #" + std::to_string(i+1) + ").");
        }
        start = end+1;
    }

    if (stop) *stop = end - 1;
    return ngpt::datetime<T> {year{ints[0]}, month{ints[1]}, day_of_month{ints[2]},
        hours{ints[3]}, minutes{ints[4]}, T{0}};
}

} // namespace igs_log_details

class igs_log
{
public:
    explicit
    igs_log(const std::string& filename)
    : m_filename{filename},
      m_ifs{filename.c_str(), std::ifstream::in}
    {
        if ( !m_ifs.is_open() ) {
            throw std::invalid_argument
            ("igs_log: Could not open log file: \""+filename+"\"");
        }
    }

    /// No copy constructor.
    igs_log(const igs_log&) = delete;

    /// No assignment operator.
    igs_log& operator=(const igs_log&) = delete;

    /// Move constructor.
    igs_log(igs_log&& ec)
    : m_filename{std::move(ec.m_filename)},
      m_ifs{std::move(ec.m_ifs)}
    {}

    /// Move assignment operator.
    igs_log& operator=(igs_log&& ec)
    {
        if (this != &ec) {
            m_filename = std::move(ec.m_filename);
            m_ifs = std::move(ec.m_ifs);
        }
        return *this;
    }

    template<typename T>
        void
        receiver_changes()
    {
        constexpr std::size_t MAX_CHARS {256};
        const char* receiver_info_block = "3.   GNSS Receiver Information";

        char line[MAX_CHARS];
        ngpt::datetime<T> from,
                          to;
                          //prev_epoch = ngpt::datetime<T>::min();
        std::string rectype;

        this->rewind();
        
        // read all lines untill we reach the block: '3.   GNSS Receiver Information'
        while ( m_ifs.getline(line, MAX_CHARS) && strcmp(receiver_info_block, line) ) {
            std::cout<<"\n[DEBUG] Skipping line:"<<line;
        }
        
        if ( strcmp(receiver_info_block, line) ) {
            throw std::invalid_argument
                ("receiver_changes: Failed to find receiver info block in log file.");
        }

        int answer;
        while ( ! (answer = this->read_receiver_block(rectype, from, to)) ) {
            std::cout<<"\n[DEBUG] Read new block:";
            std::cout<<"\nfrom "<<ngpt::strftime_ymd_hmfs(from)<<" to "<< ngpt::strftime_ymd_hmfs(to)<<", rec: "<<rectype;
        }

        if (answer < 0 ) {
            std::cout<<"\n[ERROR] Got back "<< answer;
        }
    }

private:
    /// The name of the log file.
    std::string m_filename;
    /// The input file stream
    std::ifstream m_ifs;
    
    /// @brief Go to the begining of the stream (input file).
    void rewind() noexcept
    {
        m_ifs.seekg(0, std::ios::beg);
        return;
    }

    /// First line to be read, should be: '3.x  Receiver Type            : (A20, from rcvr_ant.tab; see instructions)'
    /// < 0 error
    /// = 0 all ok
    /// > 0 block 3.x (skipped)
    template<typename T>
        int
        read_receiver_block(std::string& receiver_type, ngpt::datetime<T>& start, ngpt::datetime<T>& stop)
    {
        constexpr std::size_t MAX_CHARS {256};
        char line[MAX_CHARS];
        m_ifs.getline(line, MAX_CHARS); // empty line
        m_ifs.getline(line, MAX_CHARS); // first receiver-block line

        if ( line[0] != '3' || line[1] != '.' ) return -1;
        if ( line[2] == 'x' ) return 1;

        if ( strncmp("Receiver Type            :", line+5, 26) ) return -1;
        receiver_type.assign(line+31, 20);

        if ( !m_ifs.getline(line, MAX_CHARS) || strncmp("Satellite System         :", line+5, 26) )
            return -1;
        if ( !m_ifs.getline(line, MAX_CHARS) || strncmp("Serial Number            :", line+5, 26) )
            return -1;
        if ( !m_ifs.getline(line, MAX_CHARS) || strncmp("Firmware Version         :", line+5, 26) )
            return -1;
        if ( !m_ifs.getline(line, MAX_CHARS) || strncmp("Elevation Cutoff Setting :", line+5, 26) )
            return -1;
        if ( !m_ifs.getline(line, MAX_CHARS) || strncmp("Date Installed           :", line+5, 26) ) {
            return -1;
        } else {
            if ( igs_log_details::str_is_empty(line+31) ) {
                start = ngpt::datetime<T>::min();
            } else {
                int pos = 31;
                while (line[pos] && line[pos] == ' ') ++pos;
                std::cout<<"\n\tHere is the string: ["<<line+pos<<"]";
                if (   !strncmp(line+pos, "(CCYY-MM-DDThh:mmZ)", 19)
                    || !strncmp(line+pos, "CCYY-MM-DDThh:mmZ", 17) ) {
                    start = ngpt::datetime<T>::min();
                } else {
                    start = igs_log_details::strptime_log<T>(line+31);
                }
            }
        }
        if ( !m_ifs.getline(line, MAX_CHARS) || strncmp("Date Removed             :", line+5, 26) ) {
            return -1;
        } else {
            if ( igs_log_details::str_is_empty(line+31) ) {
                stop = ngpt::datetime<T>::max();
            } else {
                int pos = 31;
                while (line[pos] && line[pos] == ' ') ++pos;
                std::cout<<"\n\tHere is the string: ["<<line+pos<<"]";
                if (   !strncmp(line+pos, "(CCYY-MM-DDThh:mmZ)", 19)
                    || !strncmp(line+pos, "CCYY-MM-DDThh:mmZ", 17) ) {
                    stop = ngpt::datetime<T>::max();
                } else {
                    stop = igs_log_details::strptime_log<T>(line+31);
                }
            }
        }
        if ( !m_ifs.getline(line, MAX_CHARS) || strncmp("Temperature Stabiliz.    :", line+5, 26) )
            return -1;
        if ( !m_ifs.getline(line, MAX_CHARS) || strncmp("Additional Information   :", line+5, 26) )
            return -1;
        // read additional information lines
        bool new_line = true;
        auto pos      = m_ifs.tellg();
        while ( new_line && m_ifs.getline(line, MAX_CHARS) ) {
            char* c = line;
            int  sz = strlen(line);
            int   i = 0;
            while (c && (i < sz)) {
                if (c[i] != ' ') {
                    new_line = true;
                    pos =  m_ifs.tellg();
                    break;
                }
                ++i;
            }
            new_line = false;
        }
        m_ifs.seekg(pos);

        return 0;
    }

}; // class igs_log

} // namespace ngpt

#endif
