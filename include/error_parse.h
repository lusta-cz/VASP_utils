#ifndef ERROR_PARSE_H_INCLUDED
#define ERROR_PARSE_H_INCLUDED

#include <fstream>
#include <iostream>
#include <string>

/// @brief Category of a file-parsing error.
enum class ParseErrorKind { Io, Parse, Semantic };

/**
 * @brief Carries a single parse error: its category, line number, and message.
 *
 * Populated by @c fail() and consumed by @c reportParseError().
 */
struct ParseError {
    ParseErrorKind kind{ParseErrorKind::Parse};
    int line_number{0};
    std::string message;
};

/** @brief Fill @p error and return false; used for early-return error propagation. */
inline bool fail(ParseError& error, ParseErrorKind kind, int line_number, std::string message) {
    error.kind = kind;
    error.line_number = line_number;
    error.message = std::move(message);
    return false;
}

/** @brief Print a human-readable error message to stderr. */
inline void reportParseError(const std::string& filename, const ParseError& error) {
    const char* category = "Parse error";
    if (error.kind == ParseErrorKind::Io) {
        category = "I/O error";
    } else if (error.kind == ParseErrorKind::Semantic) {
        category = "Semantic error";
    }
    std::cerr << category << " while reading " << filename;
    if (error.line_number > 0) {
        std::cerr << " at line " << error.line_number;
    }
    std::cerr << ": " << error.message << "\n";
}

/**
 * @brief Read the next line from @p file, incrementing @p line_number.
 * @return false (with @p error populated) on EOF or I/O failure.
 */
inline bool readRequiredLine(std::ifstream& file, std::string& line, int& line_number, ParseError& error,
                             const std::string& context) {
    if (!std::getline(file, line)) {
        return fail(error, file.bad() ? ParseErrorKind::Io : ParseErrorKind::Parse, line_number + 1,
                    "unexpected end of file while reading " + context);
    }
    ++line_number;
    return true;
}

#endif  // ERROR_PARSE_H_INCLUDED
