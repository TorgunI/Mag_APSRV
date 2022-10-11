#include "welding.h"



void welding_engine(struct WeldingEngine* engine, enum WeldingInputSymbol symbol)
{
	switch (engine->state)
	{
		
		case WELDING_STATE_ZERO:
		{
			if (INPUT_SYMBOL_ZERO == symbol)
			{
				engine->gaz = GAZ_ZERO;
				return;
			}
			else if (INPUT_SYMBOL_ONE == symbol)
			{
				flag = false;
				engine->gaz = GAZ_HIGH;
				engine->state = WELDING_STATE_ONE;

				return;
			}
			break;
		}
		
		case WELDING_STATE_ONE:
		{
			if (INPUT_SYMBOL_ZERO == symbol)
			{
				engine->current = CURRENT_ZERO;
				engine->arcv = ARC_ZERO;
				engine->speed = SPEED_ZERO;
				engine->oscv = OSC_ZERO;
				engine->state = WELDING_STATE_ZERO;
				return;
			}
			
			else if (INPUT_SYMBOL_ONE == symbol)
			{
				if (flag)
				{
					engine->oscv = OSC_ZERO;
				}
				else
				{
					engine->oscv = OSC_HIGH;
				}
				engine->current = HIGH_TOK;
				engine->arcv = HIGH_DUGA;
				engine->speed = HIGH_SPEED;
				flag = true;
				return;
			}
			break;
	}

	return;
	}
	
}

bool welding_init(struct WeldingEngine* engine)
	{
		if (0 == engine)
		{
			return false;
		}
		engine->state = DEFAULT_WELDING_STATE;
		return true;
	}
bool welding_reset(struct WeldingEngine* engine)
	{
		return welding_init(engine);
	}


