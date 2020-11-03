#ifndef DEMOTRACK_BEAM_ELEMENTS_H__
#define DEMOTRACK_BEAM_ELEMENTS_H__

typedef enum
{
    BEAM_ELEMENT_DRIFT = 3,
    BEAM_ELEMENT_MULTIPOLE = 4,
    BEAM_ELEMENT_CAVITY = 5,
    BEAM_ELEMENT_SC_COASTING = 16
}
beam_element_type;

typedef struct Drift
{
    double type_id;
    double length;
}
Drift;

typedef struct Multipole
{
    double type_id;
    double order;
    double length;
    double bal[ 16 ];
}
Multipole;

typedef struct Cavity
{
    double type_id;
    double voltage;
    double frequency;
    double lag;
}
Cavity;

typedef struct SpaceChargeCoasting
{
    double type_id;
    double number_of_particles;
    double circumference;
    double sigma_x;
    double sigma_y;
    double length;
    double x_co;
    double y_co;
    double min_sigma_diff;
    double enabled;
}
SpaceChargeCoasting;

typedef enum
{
    NUM_SLOTS_DRIFT = 2,
    NUM_SLOTS_MULTIPOLE = 19,
    NUM_SLOTS_CAVITY = 4,
    NUM_SLOTS_SC_COASTING = 10
}
beam_element_nslots;


#endif /* DEMOTRACK_BEAM_ELEMENTS_H__ */
