#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#define min(a, b) (a > b ? b : a)
#define max(a, b) (a < b ? b : a)

const int MAP_SIZE = 21000000;
const int WIDTH = 7000;
const int HEIGHT = 3000;

const int POPULATION_SIZE = 40;
const int ACTION_DURATION = 5;
const int ELITISM = 4;

double get_rand() {
  return rand() / (double) RAND_MAX;
}

int get_random(int min, int max) {
  return (int) roundf(get_rand() * (max - min + 1)) + min;;
}

double get_millis() {
  return ((double) clock() / (double) CLOCKS_PER_SEC) * 1000.0;
}

typedef struct {
  int x;
  int y;
} Coordinate;

Coordinate zero_coordinate() {
  return (Coordinate) {0, 0};
}

typedef struct {
  int x;
  int y;
} Vector2D;

Vector2D zero_vector() {
  return (Vector2D) {0, 0};
}

bool is_valid_coordinate(int x, int y) {
  return x >= 0 && x < WIDTH && y >= 0 && y < HEIGHT;
}

typedef struct {
  bool * ground; // 7000 * 3000
  Coordinate landing[2];
} Surface;

Surface * create_surface() {
  Surface * surface = (Surface *) malloc(sizeof(Surface));
  surface->ground = (bool *) calloc(sizeof(bool), MAP_SIZE);
  surface->landing[0].x = 0;
  surface->landing[0].y = 0;
  surface->landing[1].x = 0;
  surface->landing[1].y = 0;
  return surface;
}

typedef struct {
  int thrust;
  int angle;
  double fitness;
} Action;

Action zero_action() {
  return (Action) {0, 0, -1.0};
}

typedef struct {
  Action * actions;
} Population;

typedef struct {
  Coordinate coordinate;
  Vector2D speed;
  int fuel;
  Action action;
} Lander;

Lander initialize_lander() {
  return (Lander) {zero_coordinate(), zero_vector(), 0, zero_action()};
}

Lander advance(Lander lander, Action action) {
  float xAcc = ((float) lander.action.thrust) * sinf((float) lander.action.angle);
  float yAcc = ((float) lander.action.thrust) * cosf((float) lander.action.angle) - 3.711f;

  int nextX = (int) roundf((float) lander.coordinate.x + (float) lander.speed.x + (0.5f * xAcc));
  int nextY = (int) roundf((float) lander.coordinate.y + (float) lander.speed.y + (0.5f * yAcc));

  Coordinate next_coordinate = (Coordinate) {nextX, nextY};
  Vector2D next_speed = (Vector2D) {lander.speed.x + (int) roundf(xAcc), lander.speed.y + (int) roundf(yAcc)};
  Action next_action = (Action) { lander.action.thrust + action.thrust, lander.action.angle + action.angle, -1.0 };
  next_action.thrust = max(0, min(4, next_action.thrust));
  next_action.angle = max(-90, min(90, next_action.angle));

  Lander result = (Lander) {
    next_coordinate,
    next_speed,
    lander.fuel - lander.action.thrust,
    next_action
  };

  return result;
}

void push(int * stack, int * top, int x, int y) {
  stack[++(*top)] = x;
  stack[++(*top)] = y;
}

bool pop(const int * stack, int * top, int * x, int * y) {
  if (*top < 0) return false;
  *y = stack[(*top)--];
  *x = stack[(*top)--];
  return true;
}

void draw_line(Surface * surface, Coordinate from, Coordinate to) {
  int dx = to.x - from.x;
  int dy = to.y - from.y;
  int steps = abs(dx) > abs(dy) ? abs(dx) : abs(dy);

  float xIncrement = (float) dx / (float) steps;
  float yIncrement = (float) dy / (float) steps;

  float x = (float) from.x;
  float y = (float) from.y;

  for (int i = 0; i <= steps; i++) {
    int index = (int) (roundf(y) * (float) WIDTH + roundf(x));
    surface->ground[index] = true;
    x += xIncrement;
    y += yIncrement;
  }
}

void add_point(Surface * surface, Coordinate * coordinates, int current_position) {
  Coordinate next_point = coordinates[current_position];
  int index = next_point.y * WIDTH + next_point.x;
  surface->ground[index] = true;
  if (current_position > 0) {
    Coordinate last_point = coordinates[current_position - 1];
    if (last_point.y == next_point.y) {
      surface->landing[0].x = last_point.x;
      surface->landing[0].y = last_point.y;
      surface->landing[1].x = next_point.x;
      surface->landing[1].y = next_point.y;
    }
    draw_line(surface, last_point, next_point);
  }
}

void flood_fill_scanline(Surface * surface, Coordinate from) {
  int x1;

  int x = from.x;
  int y = from.y;

  int spanAbove, spanBelow;
  int * stack = malloc(42000000 * sizeof(int));

  int top = -1;

  push(stack, &top, x, y);

  while (pop(stack, &top, &x, &y)) {
    x1 = x;
    while (x1 >= 0 && is_valid_coordinate(x1, y) && !surface->ground[y * WIDTH + x1]) x1--;
    x1 += 1;
    spanAbove = spanBelow = 0;

    while (x1 < WIDTH && !surface->ground[y * WIDTH + x1]) {
      surface->ground[y * WIDTH + x1] = true;
      if(!spanAbove && y > 0 && !surface->ground[(y - 1) * WIDTH + x1]) {
        push(stack, &top, x1, y - 1);
        spanAbove = 1;
      }
      else if(spanAbove && y > 0 && surface->ground[(y - 1) * WIDTH + x1]) spanAbove = 0;
      if(!spanBelow && y < HEIGHT - 1 && !surface->ground[(y + 1) * WIDTH + x1]) {
        push(stack, &top, x1, y + 1);
        spanBelow = 1;
      }
      else if(spanBelow && y < HEIGHT - 1 && surface->ground[(y + 1) * WIDTH + x1]) spanBelow = 0;
      x1 += 1;
    }
  }

  free(stack);
}

void fill_surface(Surface * surface) {
  flood_fill_scanline(surface, (Coordinate) { (surface->landing[0].x + surface->landing[1].x) / 2, surface->landing[1].y - 5});
}

bool is_crash(Coordinate coordinate, Surface * surface) {
  int index = coordinate.y * WIDTH + coordinate.x;
  return surface->ground[index];
}

bool is_lost(Coordinate coordinate) {
  return coordinate.x < 0 || coordinate.x >= WIDTH || coordinate.y < 0 || coordinate.y >= HEIGHT;
}

bool is_game_over(Coordinate coordinate, Surface * surface) {
  return is_lost(coordinate) || is_crash(coordinate, surface);
}

double evaluate_action(Action action, Lander lander, Surface * surface) {
  Lander evaluatedLander = lander;

  int duration = 0;
  double score = -10.0;
  while (duration < ACTION_DURATION || is_game_over(evaluatedLander.coordinate, surface)) {
    score += 10;
    evaluatedLander = advance(evaluatedLander, action);
    duration += 1;
    if (is_game_over(lander.coordinate, surface)) break;
  }

//  if (grid.collidesAt(state.coordinates)) {
//                    score = if (grid.isLandingZone(state.coordinates)) {
//                        if (debug) {
//                            System.err.println("collision at $state")
//                            System.err.println("fuel: ${2 * state.fuel}")
//                            System.err.println("vSpeedPenalty: ${-max(0.0, abs(state.vSpeed) - safeVerticalSpeed) * 200}")
//                            System.err.println("hSpeedPenalty: ${-max(0.0, abs(state.hSpeed) - safeHorizontalSpeed) * 200}")
//                            System.err.println("rotatePenalty: ${-abs(state.command.rotate) * 1000}")
//                            System.err.println("targetPenalty: ${-state.coordinates.distanceFrom(grid.landingTarget) * 5}")
//                        }
//                        100000.0 + 2 * state.fuel -
//                            max(0.0, abs(state.vSpeed) - safeVerticalSpeed) * 200 -
//                            max(0.0, abs(state.hSpeed) - safeHorizontalSpeed) * 200 -
//                            abs(state.command.rotate) * 1000 -
//                            state.coordinates.distanceFrom(grid.landingTarget) * 5
//                    } else {
//                        -min(state.coordinates.distanceFrom(grid.landing.first), state.coordinates.distanceFrom(grid.landing.second))
//                    }
//                    return
//                }

//  score +=

  return score;
}

Population * initialize_population(Lander lander, Surface * surface) {
  Population * population = (Population *) malloc(sizeof(Population));
  population->actions = (Action *) calloc(93, sizeof(Action));

  int next = 0;
  for (int angle = -15; angle <= 15; angle++) {
    for (int thrust = -1; thrust <= 1; thrust++) {
      Action action = (Action) {thrust, angle};
      action.fitness = evaluate_action(action, lander, surface);
      population->actions[next] = action;
      next += 1;
    }
  }

  return population;
}

Action find_better(Population * population) {
  Action current = population->actions[0];

  for (int i = 0; i < 93; i++)
    if (population->actions[i].fitness > current.fitness)
      current = population->actions[i];

  return current;
}

int main() {
  srand(time(NULL));
  Surface * surface = create_surface();
  int surface_n;
  scanf("%d", &surface_n);

  Coordinate coordinates[surface_n + 2];
  coordinates[0] = zero_coordinate();
  add_point(surface, coordinates, 0);
  for (int i = 1; i <= surface_n; i++) {
    int land_x, land_y;
    scanf("%d%d", &land_x, &land_y);
    Coordinate coordinate = (Coordinate) {land_x, land_y};
    coordinates[i] = coordinate;
    add_point(surface, coordinates, i);
  }
  coordinates[surface_n + 1] = (Coordinate) {6999, 0};
  add_point(surface, coordinates, surface_n + 1);
  draw_line(surface, zero_coordinate(), coordinates[surface_n + 1]);

  clock_t start_time, end_time;
  start_time = clock();
  fill_surface(surface);
  end_time = clock();

  fprintf(stderr, "Filled surface in %f ms\n", ((double) (end_time - start_time) / CLOCKS_PER_SEC) * 1000.0);

  // game loop
  while (true) {
    Lander lander = initialize_lander();
    scanf(
        "%d%d%d%d%d%d%d",
        &lander.coordinate.x,
        &lander.coordinate.y,
        &lander.speed.x,
        &lander.speed.y,
        &lander.fuel,
        &lander.action.angle,
        &lander.action.thrust
    );

    double start = get_millis();
    Population * population = initialize_population(lander, surface);

    Action better = find_better(population);
//      qsort(population.actions, 93, sizeof(Action), )

    fprintf(stderr, "Selected in %f ms\n", get_millis() - start);

//    for (int i = 0; i < 93; i++) {
//      Action action = population->actions[i];
//      fprintf(stderr, "{%f}", action.fitness);
//    }

    better.thrust += lander.action.thrust;
    better.angle += lander.action.angle;
    better.thrust = min(4, max(0, better.thrust));
    better.angle = min(90, max(-90, better.angle));
    free(population->actions);
    free(population);
    printf("%d %d\n", better.angle, better.thrust);
  }

  return 0;
}
