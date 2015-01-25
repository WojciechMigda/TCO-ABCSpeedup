/*******************************************************************************
 *
 * Filename: json.hpp
 *
 * Description:
 *      description
 *
 *******************************************************************************
 * History:
 * --------
 * Date         Who  Ticket     Description
 * ----------   ---  ---------  ------------------------------------------------
 * 2015-01-16   xx              Initial version
 *
 ******************************************************************************/
#ifndef JSON_HPP_
#define JSON_HPP_

#include <string>
#include <map>
#include <limits>
#include <istream>
#include <vector>
#include <cassert>
#include <sstream>
#include <iomanip>

// jsonxx versioning: major.minor-extra where
// major = { number }
// minor = { number }
// extra = { 'a':alpha, 'b':beta, 'rc': release candidate, 'r': release, 's':stable }
#define JSONXX_MAJOR    "0"
#define JSONXX_MINOR    "22"
#define JSONXX_EXTRA    "a"
#define JSONXX_VERSION  JSONXX_MAJOR "." JSONXX_MINOR "-" JSONXX_EXTRA
#define JSONXX_XML_TAG  "<!-- generated by jsonxx " JSONXX_VERSION " -->"

#if __cplusplus > 199711L
#define JSONXX_COMPILER_HAS_CXX11 1
#elif defined(_MSC_VER) && _MSC_VER > 1700
#define JSONXX_COMPILER_HAS_CXX11 1
#else
#define JSONXX_COMPILER_HAS_CXX11 0
#endif

#define JSONXX_ASSERT(...) do { if( jsonxx::Assertions ) \
  jsonxx::assertion(__FILE__,__LINE__,#__VA_ARGS__,bool(__VA_ARGS__)); } while(0)

namespace jsonxx {

// Settings
enum Settings {
  // constants
  Enabled = true,
  Disabled = false,
  Permissive = true,
  Strict = false,
  // values
  Parser = Permissive,  // permissive or strict parsing
  UnquotedKeys = Disabled, // support of unquoted keys
  Assertions = Enabled  // enabled or disabled assertions (these asserts work both in DEBUG and RELEASE builds)
};

enum Format {
  JSON      = 0,     // JSON output
  JSONx     = 1,     // XML output, JSONx  format. see http://goo.gl/I3cxs
};

// Types
typedef long double Number;
typedef bool Boolean;
typedef std::string String;
struct Null {};
class Value;
class Object;
class Array;

// Identity meta-function
template <typename T>
struct identity {
  typedef T type;
};

// Detail
void assertion( const char *file, int line, const char *expression, bool result );

// A JSON Object
class Object {
 public:
  Object();
  ~Object();

  template <typename T>
  bool has(const std::string& key) const;

  // Always call has<>() first. If the key doesn't exist, consider
  // the behavior undefined.
  template <typename T>
  T& get(const std::string& key);
  template <typename T>
  const T& get(const std::string& key) const;

  template <typename T>
  const T& get(const std::string& key, const typename identity<T>::type& default_value) const;

  size_t size() const;
  bool empty() const;

  const std::map<std::string, Value*>& kv_map() const;

  void reset();
  bool parse(std::istream &input);
  bool parse(const std::string &input);
  typedef std::map<std::string, Value*> container;
  void import( const Object &other );
  void import( const std::string &key, const Value &value );
  Object &operator<<(const Value &value);
  Object &operator<<(const Object &value);
  Object &operator=(const Object &value);
  Object(const Object &other);
  Object(const std::string &key, const Value &value);
  template<size_t N>
  Object(const char (&key)[N], const Value &value) {
    import(key,value);
  }
  template<typename T>
  Object &operator<<(const T &value);

 protected:
  static bool parse(std::istream& input, Object& object);
  container value_map_;
  std::string odd;
};

class Array {
 public:
  Array();
  ~Array();

  size_t size() const;
  bool empty() const;

  template <typename T>
  bool has(unsigned int i) const;

  template <typename T>
  T& get(unsigned int i);
  template <typename T>
  const T& get(unsigned int i) const;

  template <typename T>
  const T& get(unsigned int i, const typename identity<T>::type& default_value) const;

  const std::vector<Value*>& values() const {
    return values_;
  }

  void reset();
  bool parse(std::istream &input);
  bool parse(const std::string &input);
  typedef std::vector<Value*> container;
  void import(const Array &other);
  void import(const Value &value);
  Array &operator<<(const Array &other);
  Array &operator<<(const Value &value);
  Array &operator=(const Array &other);
  Array &operator=(const Value &value);
  Array(const Array &other);
  Array(const Value &value);
 protected:
  static bool parse(std::istream& input, Array& array);
  container values_;
};

// A value could be a number, an array, a string, an object, a
// boolean, or null
class Value {
 public:

  Value();
  ~Value() { reset(); }
  void reset();

  template<typename T>
  void import( const T & ) {
    reset();
    type_ = INVALID_;
    // debug
    // std::cout << "[WARN] No support for " << typeid(t).name() << std::endl;
  }
  void import( const bool &b ) {
    reset();
    type_ = BOOL_;
    bool_value_ = b;
  }
#define $number(TYPE) \
  void import( const TYPE &n ) { \
    reset(); \
    type_ = NUMBER_; \
    number_value_ = static_cast<long double>(n); \
  }
  $number( char )
  $number( int )
  $number( long )
  $number( long long )
  $number( unsigned char )
  $number( unsigned int )
  $number( unsigned long )
  $number( unsigned long long )
  $number( float )
  $number( double )
  $number( long double )
#undef $number
#if JSONXX_COMPILER_HAS_CXX11 > 0
  void import( const std::nullptr_t & ) {
    reset();
    type_ = NULL_;
  }
#endif
  void import( const Null & ) {
    reset();
    type_ = NULL_;
  }
  void import( const String &s ) {
    reset();
    type_ = STRING_;
    *( string_value_ = new String() ) = s;
  }
  void import( const Array &a ) {
    reset();
    type_ = ARRAY_;
    *( array_value_ = new Array() ) = a;
  }
  void import( const Object &o ) {
    reset();
    type_ = OBJECT_;
    *( object_value_ = new Object() ) = o;
  }
  void import( const Value &other ) {
    if (this != &other)
    switch (other.type_) {
      case NULL_:
        import( Null() );
        break;
      case BOOL_:
        import( other.bool_value_ );
        break;
      case NUMBER_:
        import( other.number_value_ );
        break;
      case STRING_:
        import( *other.string_value_ );
        break;
      case ARRAY_:
        import( *other.array_value_ );
        break;
      case OBJECT_:
        import( *other.object_value_ );
        break;
      case INVALID_:
        type_ = INVALID_;
        break;
      default:
        JSONXX_ASSERT( !"not implemented" );
    }
  }
  template<typename T>
  Value &operator <<( const T &t ) {
    import(t);
    return *this;
  }
  template<typename T>
  Value &operator =( const T &t ) {
    reset();
    import(t);
    return *this;
  }
  Value(const Value &other);
  template<typename T>
  Value( const T&t ) : type_(INVALID_) { import(t); }
  template<size_t N>
  Value( const char (&t)[N] ) : type_(INVALID_) { import( std::string(t) ); }

  bool parse(std::istream &input);
  bool parse(const std::string &input);

  template<typename T>
  bool is() const;
  template<typename T>
  T& get();
  template<typename T>
  const T& get() const;

  bool empty() const;

 public:
  enum {
    NUMBER_,
    STRING_,
    BOOL_,
    NULL_,
    ARRAY_,
    OBJECT_,
    INVALID_
  } type_;
  union {
    Number number_value_;
    String* string_value_;
    Boolean bool_value_;
    Array* array_value_;
    Object* object_value_;
  };

protected:
  static bool parse(std::istream& input, Value& value);
};

template <typename T>
bool Array::has(unsigned int i) const {
  if (i >= size()) {
    return false;
  } else {
    Value* v = values_.at(i);
    return v->is<T>();
  }
}

template <typename T>
T& Array::get(unsigned int i) {
  JSONXX_ASSERT(i < size());
  Value* v = values_.at(i);
  return v->get<T>();
}

template <typename T>
const T& Array::get(unsigned int i) const {
  JSONXX_ASSERT(i < size());
  const Value* v = values_.at(i);
  return v->get<T>();
}

template <typename T>
const T& Array::get(unsigned int i, const typename identity<T>::type& default_value) const {
  if(has<T>(i)) {
    const Value* v = values_.at(i);
    return v->get<T>();
  } else {
    return default_value;
  }
}

template <typename T>
bool Object::has(const std::string& key) const {
  container::const_iterator it(value_map_.find(key));
  return it != value_map_.end() && it->second->is<T>();
}

template <typename T>
T& Object::get(const std::string& key) {
  JSONXX_ASSERT(has<T>(key));
  return value_map_.find(key)->second->get<T>();
}

template <typename T>
const T& Object::get(const std::string& key) const {
  JSONXX_ASSERT(has<T>(key));
  return value_map_.find(key)->second->get<T>();
}

template <typename T>
const T& Object::get(const std::string& key, const typename identity<T>::type& default_value) const {
  if (has<T>(key)) {
    return value_map_.find(key)->second->get<T>();
  } else {
    return default_value;
  }
}

template<>
inline bool Value::is<Value>() const {
    return true;
}

template<>
inline bool Value::is<Null>() const {
  return type_ == NULL_;
}

template<>
inline bool Value::is<Boolean>() const {
  return type_ == BOOL_;
}

template<>
inline bool Value::is<String>() const {
  return type_ == STRING_;
}

template<>
inline bool Value::is<Number>() const {
  return type_ == NUMBER_;
}

template<>
inline bool Value::is<Array>() const {
  return type_ == ARRAY_;
}

template<>
inline bool Value::is<Object>() const {
  return type_ == OBJECT_;
}

template<>
inline Value& Value::get<Value>() {
    return *this;
}

template<>
inline const Value& Value::get<Value>() const {
    return *this;
}

template<>
inline bool& Value::get<Boolean>() {
  JSONXX_ASSERT(is<Boolean>());
  return bool_value_;
}

template<>
inline std::string& Value::get<String>() {
  JSONXX_ASSERT(is<String>());
  return *string_value_;
}

template<>
inline Number& Value::get<Number>() {
  JSONXX_ASSERT(is<Number>());
  return number_value_;
}

template<>
inline Array& Value::get<Array>() {
  JSONXX_ASSERT(is<Array>());
  return *array_value_;
}

template<>
inline Object& Value::get<Object>() {
  JSONXX_ASSERT(is<Object>());
  return *object_value_;
}

template<>
inline const Boolean& Value::get<Boolean>() const {
  JSONXX_ASSERT(is<Boolean>());
  return bool_value_;
}

template<>
inline const String& Value::get<String>() const {
  JSONXX_ASSERT(is<String>());
  return *string_value_;
}

template<>
inline const Number& Value::get<Number>() const {
  JSONXX_ASSERT(is<Number>());
  return number_value_;
}

template<>
inline const Array& Value::get<Array>() const {
  JSONXX_ASSERT(is<Array>());
  return *array_value_;
}

template<>
inline const Object& Value::get<Object>() const {
  JSONXX_ASSERT(is<Object>());
  return *object_value_;
}

template<typename T>
inline Object &Object::operator<<(const T &value) {
  return *this << Value(value), *this;
}

}  // namespace jsonxx

#if defined(NDEBUG) || defined(_NDEBUG)
#   define JSONXX_REENABLE_NDEBUG
#   undef  NDEBUG
#   undef _NDEBUG
#endif
void jsonxx::assertion( const char *file, int line, const char *expression, bool result ) {
    if( !result ) {
        fprintf( stderr, "[JSONXX] expression '%s' failed at %s:%d -> ", expression, file, line );
        assert( 0 );
    }
}
#if defined(JSONXX_REENABLE_NDEBUG)
#   define  NDEBUG
#   define _NDEBUG
#endif
namespace jsonxx {

//static_assert( sizeof(unsigned long long) < sizeof(long double), "'long double' cannot hold 64bit values in this compiler :(");

bool match(const char* pattern, std::istream& input);
bool parse_array(std::istream& input, Array& array);
bool parse_bool(std::istream& input, Boolean& value);
bool parse_comment(std::istream &input);
bool parse_null(std::istream& input);
bool parse_number(std::istream& input, Number& value);
bool parse_object(std::istream& input, Object& object);
bool parse_string(std::istream& input, String& value);
bool parse_identifier(std::istream& input, String& value);
bool parse_value(std::istream& input, Value& value);

// Try to consume characters from the input stream and match the
// pattern string.
bool match(const char* pattern, std::istream& input) {
    input >> std::ws;
    const char* cur(pattern);
    char ch(0);
    while(input && !input.eof() && *cur != 0) {
        input.get(ch);
        if (ch != *cur) {
            input.putback(ch);
            if( parse_comment(input) )
                continue;
            while (cur > pattern) {
                cur--;
                input.putback(*cur);
            }
            return false;
        } else {
            cur++;
        }
    }
    return *cur == 0;
}

bool parse_string(std::istream& input, String& value) {
    char ch = '\0', delimiter = '"';
    if (!match("\"", input))  {
        if (Parser == Strict) {
            return false;
        }
        delimiter = '\'';
        if (input.peek() != delimiter) {
            return false;
        }
        input.get(ch);
    }
    while(!input.eof() && input.good()) {
        input.get(ch);
        if (ch == delimiter) {
            break;
        }
        if (ch == '\\') {
            input.get(ch);
            switch(ch) {
                case '\\':
                case '/':
                    value.push_back(ch);
                    break;
                case 'b':
                    value.push_back('\b');
                    break;
                case 'f':
                    value.push_back('\f');
                    break;
                case 'n':
                    value.push_back('\n');
                    break;
                case 'r':
                    value.push_back('\r');
                    break;
                case 't':
                    value.push_back('\t');
                    break;
                case 'u': {
                        int i;
                        std::stringstream ss;
                        for( i = 0; (!input.eof() && input.good()) && i < 4; ++i ) {
                            input.get(ch);
                            ss << std::hex << ch;
                        }
                        if( input.good() && (ss >> i) )
                            value.push_back(i);
                    }
                    break;
                default:
                    if (ch != delimiter) {
                        value.push_back('\\');
                        value.push_back(ch);
                    } else value.push_back(ch);
                    break;
            }
        } else {
            value.push_back(ch);
        }
    }
    if (input && ch == delimiter) {
        return true;
    } else {
        return false;
    }
}

bool parse_identifier(std::istream& input, String& value) {
    input >> std::ws;

    char ch = '\0', delimiter = ':';
    bool first = true;

    while(!input.eof() && input.good()) {
        input.get(ch);

        if (ch == delimiter) {
            input.unget();
            break;
        }

        if(first) {
            if ((ch != '_' && ch != '$') &&
                    (ch < 'a' || ch > 'z') &&
                    (ch < 'A' || ch > 'Z')) {
                return false;
            }
            first = false;
        }
        if(ch == '_' || ch == '$' ||
            (ch >= 'a' && ch <= 'z') ||
            (ch >= 'A' && ch <= 'Z') ||
            (ch >= '0' && ch <= '9')) {
            value.push_back(ch);
        }
        else if(ch == '\t' || ch == ' ') {
            input >> std::ws;
        }
    }
    if (input && ch == delimiter) {
        return true;
    } else {
        return false;
    }
}

bool parse_number(std::istream& input, Number& value) {
    input >> std::ws;
    std::streampos rollback = input.tellg();
    input >> value;
    if (input.fail()) {
        input.clear();
        input.seekg(rollback);
        return false;
    }
    return true;
}

bool parse_bool(std::istream& input, Boolean& value) {
    if (match("true", input))  {
        value = true;
        return true;
    }
    if (match("false", input)) {
        value = false;
        return true;
    }
    return false;
}

bool parse_null(std::istream& input) {
    if (match("null", input))  {
        return true;
    }
    if (Parser == Strict) {
        return false;
    }
    return (input.peek()==',');
}

bool parse_array(std::istream& input, Array& array) {
    return array.parse(input);
}

bool parse_object(std::istream& input, Object& object) {
    return object.parse(input);
}

bool parse_comment(std::istream &input) {
    if( Parser == Permissive )
    if( !input.eof() && input.peek() == '/' )
    {
        char ch0(0);
        input.get(ch0);

        if( !input.eof() )
        {
            char ch1(0);
            input.get(ch1);

            if( ch0 == '/' && ch1 == '/' )
            {
                // trim chars till \r or \n
                for( char ch(0); !input.eof() && (input.peek() != '\r' && input.peek() != '\n'); )
                    input.get(ch);

                // consume spaces, tabs, \r or \n, in case no eof is found
                if( !input.eof() )
                    input >> std::ws;
                return true;
            }

            input.unget();
            input.clear();
        }

        input.unget();
        input.clear();
    }

    return false;
}

bool parse_value(std::istream& input, Value& value) {
    return value.parse(input);
}


Object::Object() : value_map_() {}

Object::~Object() {
    reset();
}

bool Object::parse(std::istream& input, Object& object) {
    object.reset();

    if (!match("{", input)) {
        return false;
    }
    if (match("}", input)) {
        return true;
    }

    do {
        std::string key;
        if(UnquotedKeys == Enabled) {
            if (!parse_identifier(input, key)) {
                if (Parser == Permissive) {
                    if (input.peek() == '}')
                        break;
                }
                return false;
            }
        }
        else {
            if (!parse_string(input, key)) {
                if (Parser == Permissive) {
                    if (input.peek() == '}')
                        break;
                }
                return false;
            }
        }
        if (!match(":", input)) {
            return false;
        }
        Value* v = new Value();
        if (!parse_value(input, *v)) {
            delete v;
            break;
        }
        object.value_map_[key] = v;
    } while (match(",", input));


    if (!match("}", input)) {
        return false;
    }

    return true;
}

Value::Value() : type_(INVALID_) {}

void Value::reset() {
    if (type_ == STRING_) {
        delete string_value_;
        string_value_ = 0;
    }
    else if (type_ == OBJECT_) {
        delete object_value_;
        object_value_ = 0;
    }
    else if (type_ == ARRAY_) {
        delete array_value_;
        array_value_ = 0;
    }
}

bool Value::parse(std::istream& input, Value& value) {
    value.reset();

    std::string string_value;
    if (parse_string(input, string_value)) {
        value.string_value_ = new std::string();
        value.string_value_->swap(string_value);
        value.type_ = STRING_;
        return true;
    }
    if (parse_number(input, value.number_value_)) {
        value.type_ = NUMBER_;
        return true;
    }

    if (parse_bool(input, value.bool_value_)) {
        value.type_ = BOOL_;
        return true;
    }
    if (parse_null(input)) {
        value.type_ = NULL_;
        return true;
    }
    if (input.peek() == '[') {
        value.array_value_ = new Array();
        if (parse_array(input, *value.array_value_)) {
            value.type_ = ARRAY_;
            return true;
        }
        delete value.array_value_;
    }
    value.object_value_ = new Object();
    if (parse_object(input, *value.object_value_)) {
        value.type_ = OBJECT_;
        return true;
    }
    delete value.object_value_;
    return false;
}

Array::Array() : values_() {}

Array::~Array() {
    reset();
}

bool Array::parse(std::istream& input, Array& array) {
    array.reset();

    if (!match("[", input)) {
        return false;
    }
    if (match("]", input)) {
        return true;
    }

    do {
        Value* v = new Value();
        if (!parse_value(input, *v)) {
            delete v;
            break;
        }
        array.values_.push_back(v);
    } while (match(",", input));

    if (!match("]", input)) {
        return false;
    }
    return true;
}

}  // namespace jsonxx


namespace jsonxx {
namespace {

typedef unsigned char byte;

//template<bool quote>
std::string escape_string( const std::string &input, const bool quote = false ) {
    static std::string map[256], *once = 0;
    if( !once ) {
        // base
        for( int i = 0; i < 256; ++i ) {
            map[ i ] = std::string() + char(i);
        }
        // non-printable
        for( int i = 0; i < 32; ++i ) {
            std::stringstream str;
            str << "\\u" << std::hex << std::setw(4) << std::setfill('0') << i;
            map[ i ] = str.str();
        }
        // exceptions
        map[ byte('"') ] = "\\\"";
        map[ byte('\\') ] = "\\\\";
        map[ byte('/') ] = "\\/";
        map[ byte('\b') ] = "\\b";
        map[ byte('\f') ] = "\\f";
        map[ byte('\n') ] = "\\n";
        map[ byte('\r') ] = "\\r";
        map[ byte('\t') ] = "\\t";

        once = map;
    }
    std::string output;
    output.reserve( input.size() * 2 + 2 ); // worst scenario
    if( quote ) output += '"';
    for( std::string::const_iterator it = input.begin(), end = input.end(); it != end; ++it )
        output += map[ byte(*it) ];
    if( quote ) output += '"';
    return output;
}


namespace json {

    std::string remove_last_comma( const std::string &_input ) {
        std::string input( _input );
        size_t size = input.size();
        if( size > 2 )
            if( input[ size - 2 ] == ',' )
                input[ size - 2 ] = ' ';
        return input;
    }

    std::string tag( unsigned format, unsigned depth, const std::string &name, const jsonxx::Value &t) {
        std::stringstream ss;
        const std::string tab(depth, '\t');

        if( !name.empty() )
            ss << tab << '\"' << escape_string( name ) << '\"' << ':' << ' ';
        else
            ss << tab;

        switch( t.type_ )
        {
            default:
            case jsonxx::Value::NULL_:
                ss << "null";
                return ss.str() + ",\n";

            case jsonxx::Value::BOOL_:
                ss << ( t.bool_value_ ? "true" : "false" );
                return ss.str() + ",\n";

            case jsonxx::Value::ARRAY_:
                ss << "[\n";
                for(Array::container::const_iterator it = t.array_value_->values().begin(),
                    end = t.array_value_->values().end(); it != end; ++it )
                  ss << tag( format, depth+1, std::string(), **it );
                return remove_last_comma( ss.str() ) + tab + "]" ",\n";

            case jsonxx::Value::STRING_:
                ss << '\"' << escape_string( *t.string_value_ ) << '\"';
                return ss.str() + ",\n";

            case jsonxx::Value::OBJECT_:
                ss << "{\n";
                for(Object::container::const_iterator it=t.object_value_->kv_map().begin(),
                    end = t.object_value_->kv_map().end(); it != end ; ++it)
                  ss << tag( format, depth+1, it->first, *it->second );
                return remove_last_comma( ss.str() ) + tab + "}" ",\n";

            case jsonxx::Value::NUMBER_:
                // max precision
                ss << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
                ss << t.number_value_;
                return ss.str() + ",\n";
        }
    }
} // namespace jsonxx::anon::json


} // namespace jsonxx::anon

Object::Object(const Object &other) {
  import(other);
}
Object::Object(const std::string &key, const Value &value) {
  import(key,value);
}
void Object::import( const Object &other ) {
  odd.clear();
  if (this != &other) {
    // default
    container::const_iterator
        it = other.value_map_.begin(),
        end = other.value_map_.end();
    for (/**/ ; it != end ; ++it) {
      container::iterator found = value_map_.find(it->first);
      if( found != value_map_.end() ) {
        delete found->second;
      }
      value_map_[ it->first ] = new Value( *it->second );
    }
  } else {
    // recursion is supported here
    import( Object(*this) );
  }
}
void Object::import( const std::string &key, const Value &value ) {
  odd.clear();
  container::iterator found = value_map_.find(key);
  if( found != value_map_.end() ) {
    delete found->second;
  }
  value_map_[ key ] = new Value( value );
}
Object &Object::operator=(const Object &other) {
  odd.clear();
  if (this != &other) {
    reset();
    import(other);
  }
  return *this;
}
Object &Object::operator<<(const Value &value) {
  if (odd.empty()) {
    odd = value.get<String>();
  } else {
    import( Object(odd, value) );
    odd.clear();
  }
  return *this;
}
Object &Object::operator<<(const Object &value) {
  import( std::string(odd),value);
  odd.clear();
  return *this;
}
size_t Object::size() const {
  return value_map_.size();
}
bool Object::empty() const {
  return value_map_.size() == 0;
}
const std::map<std::string, Value*> &Object::kv_map() const {
  return value_map_;
}
void Object::reset() {
  container::iterator i;
  for (i = value_map_.begin(); i != value_map_.end(); ++i) {
    delete i->second;
  }
  value_map_.clear();
}
bool Object::parse(std::istream &input) {
  return parse(input,*this);
}
bool Object::parse(const std::string &input) {
  std::istringstream is( input );
  return parse(is,*this);
}


Array::Array(const Array &other) {
  import(other);
}
Array::Array(const Value &value) {
  import(value);
}
void Array::import(const Array &other) {
  if (this != &other) {
    // default
    container::const_iterator
        it = other.values_.begin(),
        end = other.values_.end();
    for (/**/ ; it != end; ++it) {
      values_.push_back( new Value(**it) );
    }
  } else {
    // recursion is supported here
    import( Array(*this) );
  }
}
void Array::import(const Value &value) {
  values_.push_back( new Value(value) );
}
size_t Array::size() const {
  return values_.size();
}
bool Array::empty() const {
  return values_.size() == 0;
}
void Array::reset() {
  for (container::iterator i = values_.begin(); i != values_.end(); ++i) {
    delete *i;
  }
  values_.clear();
}
bool Array::parse(std::istream &input) {
  return parse(input,*this);
}
bool Array::parse(const std::string &input) {
  std::istringstream is(input);
  return parse(is,*this);
}
Array &Array::operator<<(const Array &other) {
  import(other);
  return *this;
}
Array &Array::operator<<(const Value &value) {
  import(value);
  return *this;
}
Array &Array::operator=(const Array &other) {
  if( this != &other ) {
    reset();
    import(other);
  }
  return *this;
}
Array &Array::operator=(const Value &value) {
  reset();
  import(value);
  return *this;
}

Value::Value(const Value &other) : type_(INVALID_) {
  import( other );
}
bool Value::empty() const {
  if( type_ == INVALID_ ) return true;
  if( type_ == STRING_ && string_value_ == 0 ) return true;
  if( type_ == ARRAY_ && array_value_ == 0 ) return true;
  if( type_ == OBJECT_ && object_value_ == 0 ) return true;
  return false;
}
bool Value::parse(std::istream &input) {
  return parse(input,*this);
}
bool Value::parse(const std::string &input) {
  std::istringstream is( input );
  return parse(is,*this);
}

}  // namespace jsonxx




#endif /* JSON_HPP_ */
