#include <string>

class igs_log
{
public:
    explicit
    igs_log(const std::string& filename)
    : m_filename{filename}
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
        if (*this != ec) {
            m_filename = std::move(ec.m_filename);
            m_ifs = std::move(ec.m_ifs);
        }
        return *this;
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
};

template<typename T>
    ngpt::datetime<T> strptime_log(const char* str)
{
    char *end;
    const char* start = str;
    int ints[5];

    for (int i = 0; i < 5; ++i) {
        ints[i] = static_cast<int>( std::abs(std::strtol(start, &end, 10)) );
        if (errno == ERANGE || start == end) {
            errno = 0;
            throw std::invalid_argument
                ("Invalid date format: \""+std::string(str)+"\" (argument #" + std::to_string(i+1) + ").");
        }
        start = end+1;
    }
}

/// < 0 error
/// = 0 all ok
/// > 0 block 3.x (skipped)
int
read_receiver_block(std::string& receiver_type,)
{
    constexpr std::size_t MAX_CHARS {256};
    char line[MAX_CHARS];
    m_ifs.getline(line, MAX_CHARS);

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
    }
    if ( !m_ifs.getline(line, MAX_CHARS) || strncmp("Date Removed             :", line+5, 26) ) {
        return -1;
    } else {
    }
    if ( !m_ifs.getline(line, MAX_CHARS) || strncmp("Temperature Stabiliz.    :", line+5, 26) )
        return -1;
    if ( !m_ifs.getline(line, MAX_CHARS) || strncmp("Additional Information   :", line+5, 26) )
        return -1;

    return 0;
}

void
    receiver_changes()
{
    constexpr std::size_t MAX_CHARS {256};
    const char* receiver_info_block = "3.   GNSS Receiver Information";

    char line[MAX_CHARS];
    char *cptr(&line[0]),
         *end;
    ngpt::datetime<T> epoch,
                      prev_epoch = ngpt::datetime<T>::min();

    this->rewind();
    // read all lines untill we reach the block: '3.   GNSS Receiver Information'
    while ( m_ifs.getline(line, MAX_CHARS) && strcmp(receiver_info_block, line) ) {
    }
    if ( strcmp(receiver_info_block, line) ) {
        throw std::invalid_argument
            ("receiver_changes: Failed to find receiver info block in log file.");
    }

}
