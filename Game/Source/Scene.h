#ifndef __SCENE_H__
#define __SCENE_H__

#include "Module.h"
#include <vector>

struct SDL_Texture;
struct connection
{
	int a, b;
};
class Scene : public Module
{
public:

	Scene();

	// Destructor
	virtual ~Scene();

	// Called before render is available
	bool Awake();

	// Called before the first frame
	bool Start();

	// Called before all Updates
	bool PreUpdate();

	// Called each loop iteration
	bool Update(float dt);

	// Called before all Updates
	bool PostUpdate();

	// Called before quitting
	bool CleanUp();

private:
	std::vector<Vector3> points{
		{100 * 2,100 * 2,100 * 2},
		{200 * 2,100 * 2,100 * 2},
		{200 * 2,200 * 2,100 * 2},
		{100 * 2,200 * 2,100 * 2},

		{100 * 2,100 * 2,200 * 2},
		{200 * 2,100 * 2,200 * 2},
		{200 * 2,200 * 2,200 * 2},
		{100 * 2,200 * 2,200 * 2}
	};
	std::vector<connection> connections
	{
		{0,4},
		{1,5},
		{2,6},
		{3,7},

		{0,1},
		{1,2},
		{2,3},
		{3,0},

		{4,5},
		{5,6},
		{6,7},
		{7,4}
	};
	Vector3 axis;
	Vector3 c;
	
};

#endif // __SCENE_H__