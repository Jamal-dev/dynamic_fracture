#ifndef TEST_CASES_H
#define TEST_CASES_H

#include <iostream>
#include <set>
#include <algorithm>


using namespace std;
class test_cases
{
public:
  enum Value : uint8_t
  {
    MIEHE_TENSION
  , MIEHE_SHEAR
  , L_SHAPED
  , SNEDDON
  ,PRESSURIZED
  ,SCREW_DOMI
  ,SNEDDON3D
  ,HET3D
  ,DYNAMIC_SLIT
  ,P_MESH_1
	, P_NOTCHED_CAVITY
	, P_ASYMMETRY
  };
  
  std::set<test_cases> all_cases()  { 
      return {MIEHE_TENSION
            ,MIEHE_SHEAR
          , L_SHAPED
          , SNEDDON
          ,PRESSURIZED
          ,SCREW_DOMI
          ,SNEDDON3D
          ,HET3D
          ,DYNAMIC_SLIT
					, P_NOTCHED_CAVITY
					, P_ASYMMETRY
          ,P_MESH_1}; 
 
    }
    
    constexpr bool Is_miehe_tension() const { return value == MIEHE_TENSION; }
    constexpr bool Is_miehe_shear() const { return value == MIEHE_SHEAR; }
    constexpr bool Is_l_shaped() const { return value == L_SHAPED; }
    constexpr bool Is_snedon() const { return value == SNEDDON; }
    constexpr bool Is_snedon3d() const { return value == SNEDDON3D; }
    constexpr bool Is_pressurized() const { return value == PRESSURIZED; }
    constexpr bool Is_screwdomi() const { return value == SCREW_DOMI; }
    constexpr bool Is_het3d() const { return value == HET3D; }
    constexpr bool Is_dynamicslit() const { return value == DYNAMIC_SLIT; }
    constexpr bool Is_pmesh1() const { return value == P_MESH_1; }
		constexpr bool Is_p_notched_cavity() const { return value == P_NOTCHED_CAVITY; }
    constexpr bool Is_p_asymmetry() const { return value == P_ASYMMETRY; }		
  
  
  
  test_cases() = default;
  
  constexpr test_cases(Value atest_case) : value(atest_case) { }

// #if Enable switch(test_cases) use case:
//   // Allow switch and comparisons.
//   constexpr operator Value() const { return value; }

//   // Prevent usage: if(test_case)
//   explicit operator bool() const = delete;        
// #else
// constexpr bool operator==(test_cases a) const { return value == a.value; }
//   constexpr bool operator!=(test_cases a) const { return value != a.value; }
// #endif
// Allow switch and comparisons.
  constexpr operator Value() const { return value; }

//   Prevent usage: if(test_case)
  explicit operator bool() const = delete; 

  friend std::ostream &operator<<(std::ostream &os, test_cases const &a) { 
    return os << a.value;}
    
  bool any(std::set<test_cases> const & remain_cases)
  {
      if (auto search = remain_cases.find(value); search != remain_cases.end() )
        return true;
      else
        return false;
  }
    
  std::set<test_cases> set_diff(std::set<test_cases> const & selected_test_cases){
      // set_diff function
      std::set<test_cases> _all_cases = this-> all_cases();
      std::set<test_cases> diff;
      std::set_difference(_all_cases.begin(), _all_cases.end(), selected_test_cases.begin(), selected_test_cases.end(),
                        std::inserter(diff, diff.begin()));
      return diff;
  }
  

  

private:
  Value value;

};
#endif
