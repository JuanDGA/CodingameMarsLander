#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>

const int MAP_SIZE = 21000000;
const int WIDTH = 7000;
const int HEIGHT = 3000;

typedef struct {
  int x;
  int y;
} Coordinate;

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
    int index = ((int) y * WIDTH) + (int) x;
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
    x1++;
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
      x1++;
    }
  }

  free(stack);
}

void fill_surface(Surface * surface) {
  flood_fill_scanline(surface, (Coordinate) { (surface->landing[0].x + surface->landing[1].x) / 2, surface->landing[1].y - 5});
}


int main() {
  Surface * surface = create_surface();
  fprintf(stderr, "Done");
  int surface_n;
  scanf("%d", &surface_n);

  Coordinate coordinates[surface_n];

  for (int i = 0; i < surface_n; i++) {
    int land_x, land_y;
    scanf("%d%d", &land_x, &land_y);
    Coordinate coordinate = (Coordinate) {land_x, land_y};
    coordinates[i] = coordinate;
    add_point(surface, coordinates, i);
  }

  clock_t start_time, end_time;
  start_time = clock();
  fill_surface(surface);
  end_time = clock();

  fprintf(stderr, "Filled surface in %f ms\n", ((double)(end_time - start_time) / CLOCKS_PER_SEC) * 1000.0);


  // game loop
  while (1) {
    int X;
    int Y;
    // the horizontal speed (in m/s), can be negative.
    int h_speed;
    // the vertical speed (in m/s), can be negative.
    int v_speed;
    // the quantity of remaining fuel in liters.
    int fuel;
    // the rotation angle in degrees (-90 to 90).
    int rotate;
    // the thrust power (0 to 4).
    int power;
    scanf("%d%d%d%d%d%d%d", &X, &Y, &h_speed, &v_speed, &fuel, &rotate, &power);

    // Write an action using printf(). DON'T FORGET THE TRAILING \n
    // To debug: fprintf(stderr, "Debug messages...\n");


    // rotate power. rotate is the desired rotation angle. power is the desired thrust power.
    printf("-20 3\n");
  }

  return 0;
}