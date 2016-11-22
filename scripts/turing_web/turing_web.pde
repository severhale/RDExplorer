class Cell {
  private boolean activated;
  private boolean hidden;
  private boolean aboveCell;
  private float concentration;
  public float r1, r2, c;
  public Cell(float r1, float r2, float c, boolean active) {
    activated = active;
    hidden = false;
    aboveCell = false;
    concentration = 0;
    this.r1 = r1;
    this.r2 = r2;
    this.c = c;
  }
  public Cell(Cell c) {
    activated = c.activated;
    hidden = false;
    aboveCell = false;
    concentration = c.concentration;
    r1 = c.r1;
    r2 = c.r2;
    this.c = c.c;
  }
  public void activate() {
    activated = true;
  }
  public void deactivate() {
    activated = false;
  }
  public boolean isActivated() {
    return activated;
  }
  public void addConcentration(float c) {
    concentration += c;
  }
  public void step() {
    activated = concentration > 0;
    concentration = 0;
  }
  public boolean isHidden() {
    return hidden;
  }
  public boolean aboveCell() {
    return aboveCell;
  }
  public void markAboveCell() {
    aboveCell = true;
  }
  public void markHidden() {
    hidden = true;
  }
  public void setR1(float val) {
    r1 = val;
  }
  public void setR2(float val) {
    r2 = val;
  }
}

class Grid {
  Cell[][] grid;
  Cell[][] newGrid;
  int w, h;
  float r1g, r2g, c;
  public Grid(int w, int h, float r1, float r2, float c) {
    this.w = w;
    this.h = h;
    this.r1g = r1;
    this.r2g = r2;
    this.c = c;
    grid = new Cell[w][h];
    newGrid = new Cell[w][h];
    for (int i=0; i<w; i++) {
      for (int j=0; j<h; j++) {
        //grid[i][j] = new Cell(lerp(r1-2, r1+2, (float)i/w), lerp(r2-2, r2+2, (float)j/h), c);
        grid[i][j] = new Cell(r1, r2, c, (abs(float(i)/w-.5)<.1 && abs(float(j)/h-.5)<.1 && random(1)<.5));
        newGrid[i][j] = new Cell(grid[i][j]);
      }
    }
  }
  public void distribute(Cell[][] target, float conc, int x, int y, float r) {
    for (int i0=int(x-r); i0<=x+r; i0++) {
      for (int j0=int(y-r); j0<=y+r; j0++) {
        float d = sqrt((x-i0)*(x-i0)+(y-j0)*(y-j0));
        if (d <= r) {
          if (i0 < 0 || i0 >= w || j0 < 0 || j0 >= h) {
          }
          else {
            target[i0][j0].addConcentration(d/r*conc);
          }
        }
      }
    }
  }
  public Cell getCell(int i, int j) {
    if (i >= 0 && i < w && j >= 0 && j < h) {
      return grid[i][j];
    }
    else {
      return new Cell(r1g, r2g, c, false);
    }
  }
  public boolean activatedOnGrid(int i, int j) {
    if (i >= 0 && i < w && j >= 0 && j < h) {
      return grid[i][j].isActivated();
    }
    else {
      return false;
    }
  }
  public void updateNewGrid() {
    for (int i=0; i<w; i++) {
      for (int j=0; j<h; j++) {
        if (activatedOnGrid(i, j)) {
          distribute(newGrid, 1, i, j, grid[i][j].r1);
          distribute(newGrid, grid[i][j].c, i, j, grid[i][j].r2);
        }
      }
    }
  }
  public void markHidden() {
    for (int i=0; i<w; i++) {
      for (int j=0; j<h; j++) {
        newGrid[i][j].step();
        if (newGrid[i][j].isActivated() && grid[i][j].isActivated()) {
          newGrid[i][j].markAboveCell();
          if (activatedOnGrid(i-1, j) && activatedOnGrid(i+1, j) && activatedOnGrid(i, j-1) && activatedOnGrid(i, j+1)) {          
            grid[i][j].markHidden();
          }
        }
      }
    }
  }
  public void switchToNewGrid() {
    for (int i=0; i<w; i++) {
      for (int j=0; j<h; j++) {
        grid[i][j] = new Cell(newGrid[i][j]);
      }
    }
  }
  
  public void draw(float thickness) {
    noFill();
    for (float i=0; i<w; i++) {
      for (float j=0; j<h; j++) {
        if (activatedOnGrid(int(i), int(j))) {
          stroke(0);
        }
        else {
          noStroke();
        }
        pushMatrix();
        translate(i/w*width, j/h*height);
        box(width/w, height/h, thickness);
        popMatrix();
      }
    }
  }
  public void setNewRadii(float r01, float r02) {
    for (int i=0; i<w; i++) {
      for (int j=0; j<h; j++) {
        Cell cel = grid[i][j];
        cel.setR1(r01);
        cel.setR2(r02);
      }
    }
  }
}

class GridBuilder {
  private Grid g;
  float z;
  float dz;
  ArrayList<GridLayer> layers;
  
  class GridLayer {
    ArrayList<PVector> boxes;
    int w, h;
    public GridLayer(Grid mygrid, float z) {
      boxes = new ArrayList<PVector>();
      w = mygrid.w;
      h = mygrid.h;
      for (int i=0; i<g.w; i++) {
        for (int j=0; j<g.h; j++) {
          if (mygrid.activatedOnGrid(i,j) && !mygrid.getCell(i, j).isHidden()) {
            boxes.add(new PVector((float)i/mygrid.w*width, (float)j/mygrid.h*height, z));
          }
        }
      }
    }
    public void draw(float thickness) {
      for (PVector b : boxes) {
        pushMatrix();
        translate(b.x, b.y, b.z);
        box(width/w, height/h, thickness);
        popMatrix();
      }
    }
    public int countBoxes() {
      return boxes.size();
    }
  }
  
  public GridBuilder(float z0, float dz, int w, int h, float r11, float r22, float c) {
    g = new Grid(w, h, r11, r22, c);
    this.z = z0;
    this.dz = dz;
    layers = new ArrayList<GridLayer>();
  }
  public void update() {
    // println("Starting update");
    g.updateNewGrid();
    // println("Updated new grid");
    if (layers.size() > 0) {
      g.markHidden();
      layers.set(layers.size()-1, new GridLayer(g, z - dz));
    }
    // println("Hid hidden blocks");
    g.switchToNewGrid();
    // println("Switched to new grid");
    layers.add(new GridLayer(g, z));
    // println("Added new layer to layers");
    z += dz;
    // println("Incremented z");
  }
  public void draw(float thickness) {
    for (GridLayer l : layers) {
      l.draw(thickness);
    }
  }
  
  public Grid getGrid() {
    return g;
  }
  public void reset(float z0, float are1, float are2) {
    z = z0;
    g = new Grid(g.w, g.h, are1, are2, g.c);
    layers = new ArrayList<GridLayer>();
  }
}

GridBuilder gb;
String folder;
float activatorRad = 4;
float inhibitorRad = 8;
float c = -.25;

void setup() {
  size(450, 450, P3D);
  folder = month() + "-" + day() + "-" + year() + "/" + hour() + ";" + minute() + "," + second();
  
  noStroke();
  fill(128);
  gb = new GridBuilder(-500, 5, 100, 100, activatorRad, inhibitorRad, c);
  background(255);
}

void draw() {
  lights();
  
  pushMatrix();
  
  background(240);
  
  translate(width/2, height/2, gb.z);
  rotateX(PI/2.5);
  rotateZ(frameCount * .01);
  rotateY(PI);
  translate(-width/2, -height/2, -(-500 + gb.z)/2);
  
  gb.draw(5);
  
  popMatrix();
}

void mouseClicked() {
  gb.update();
}

void keyPressed() {
  if (key == ' ') {
    // println("Resetting");
    //gb = new GridBuilder(-500, 5, 100, 100, activatorRad, inhibitorRad, c);
    gb.reset(-500, activatorRad, inhibitorRad);
    // println("Finished re-initializing");
  }
}

void setActivatorRadius(float r) {
  activatorRad = r;
  gb.getGrid().setNewRadii(activatorRad, inhibitorRad);
  // println("Set new activator radius to " + activatorRad);
}

void setInhibitorRadius(float r) {
  inhibitorRad = r;
  gb.getGrid().setNewRadii(activatorRad, inhibitorRad);
  // println("Set new inhibitor radius to " + inhibitorRad);
}