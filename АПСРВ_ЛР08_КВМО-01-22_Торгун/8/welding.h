#ifndef welding
#define welding

#include <stdbool.h>

#if defined(__cplusplus)
extern "C" {
#endif

enum WeldingState
{
	WELDING_STATE_ZERO,
	WELDING_STATE_ONE
};

enum WeldingInputSymbol
{
	INPUT_SYMBOL_ZERO = 0,
	INPUT_SYMBOL_ONE = 1
};

enum GazConsume
{
	GAZ_ZERO = 0,
	GAZ_HIGH = 25
};

enum WeldingCurrent
{
	CURRENT_ZERO = 0,
	CURRENT_HIGH = 140
};

enum ArcVoltage
{
	ARC_ZERO = 0,
	ARC_HIGH = 20
};

enum WeldingSpeed
{
	SPEED_ZERO = 0,
	SPEED_HIGH = 15
};

enum OscVoltage
{
	OSC_ZERO = 0,
	OSC_HIGH = 220
};

struct WeldingEngine
{
	enum WeldingState state;
	enum GazConsume gaz;
	enum WeldingCurrent current;
	enum ArcVoltage arcv;
	enum WeldingSpeed speed;
	enum OscVoltage oscv;
};
#define DEFAULT_WELDING_STATE WELDING_STATE_ZERO

static bool flag = false;


bool welding_init(struct WeldingEngine* engine);
bool welding_reset(struct WeldingEngine* engine);
void welding_engine(struct WeldingEngine* engine, enum WeldingInputSymbol symbol);

#endif


