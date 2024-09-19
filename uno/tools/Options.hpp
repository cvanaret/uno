// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIONS_H
#define UNO_OPTIONS_H

#include <map>
#include <string>

namespace uno {
   class Options {
   public:
      Options() = default;

      [[nodiscard]] std::string& operator[](const std::string& key);
      [[nodiscard]] const std::string& at(const std::string& key) const;

      [[nodiscard]] const std::string& get_string(const std::string& key) const;
      [[nodiscard]] double get_double(const std::string& key) const;
      [[nodiscard]] int get_int(const std::string& key) const;
      [[nodiscard]] size_t get_unsigned_int(const std::string& key) const;
      [[nodiscard]] bool get_bool(const std::string& key) const;

      void get_command_line_arguments(int argc, char* argv[]);

      std::string to_string(bool only_used) const;
      void print(bool only_used) const;

      static Options get_default_options(const std::string& file_name);

   private:
      std::map<std::string, std::string> options{};
      // keep track of the options that are used
      mutable std::map<std::string, bool> is_used{};
   };
} // namespace

#endif // UNO_OPTIONS_H
