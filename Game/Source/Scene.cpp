#include "App.h"
#include "Input.h"
#include "Textures.h"
#include "Render.h"
#include "Window.h"
#include "Scene.h"

#include "Defs.h"
#include "Log.h"

Scene::Scene() : Module()
{
	name.Create("scene");
}

// Destructor
Scene::~Scene()
{}

// Called before render is available
bool Scene::Awake()
{
	LOG("Loading Scene");
	bool ret = true;

	return ret;
}

// Called before the first frame
bool Scene::Start()
{
	axis.x = 0.1;
	axis.y = 0.1;
	axis.z = 0.1;

	for (auto& p : points)
	{
		c.x += p.x;
		c.y += p.y;
		c.z += p.z;
	}
	c.x /= points.size();
	c.y /= points.size();
	c.z /= points.size();
	return true;
}

// Called each loop iteration
bool Scene::PreUpdate()
{
	return true;
}

// Called each loop iteration
bool Scene::Update(float dt)
{
	if(app->input->GetKey(SDL_SCANCODE_UP) == KEY_REPEAT)
		app->render->camera.y -= 1;

	if(app->input->GetKey(SDL_SCANCODE_DOWN) == KEY_REPEAT)
		app->render->camera.y += 1;

	if(app->input->GetKey(SDL_SCANCODE_LEFT) == KEY_REPEAT)
		app->render->camera.x -= 1;

	if(app->input->GetKey(SDL_SCANCODE_RIGHT) == KEY_REPEAT)
		app->render->camera.x += 1;

	//app->render->DrawTexture(img, 380, 100);

	for (int i = 0; i < 8; i++)
	{
		points.at(i).x -= c.x;
		points.at(i).y -= c.y;
		points.at(i).z -= c.z;
		points.at(i).Rotate(0.002, 0.001, 0.004);
		points.at(i).x += c.x;
		points.at(i).y += c.y;
		points.at(i).z += c.z;
		app->render->Pixel(points.at(i).x, points.at(i).y);
	}
	for (int i = 0; i < 12; i++)
	{
		app->render->DrawLine(points[connections.at(i).a].x,
							  points[connections.at(i).a].y,
							  points[connections.at(i).b].x,
							  points[connections.at(i).b].y, 255, 255, 255);
	}

	app->render->PrintPoints();
	app->render->clear();
	SDL_Delay(30);

	return true;
}

// Called each loop iteration
bool Scene::PostUpdate()
{
	bool ret = true;

	if(app->input->GetKey(SDL_SCANCODE_ESCAPE) == KEY_DOWN)
		ret = false;

	return ret;
}

// Called before quitting
bool Scene::CleanUp()
{
	LOG("Freeing scene");

	return true;
}
