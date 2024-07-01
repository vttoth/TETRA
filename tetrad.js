// tetra.js : Tetrahedral constellation simulation
//
// Copyright (c) 2024 Viktor T. Toth
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// This version dated 2026/06/30.

/* -- helper functions if needed for debugging

console.logV = function(V)
{
  console.log({x:V.x.toString(),y:V.y.toString(),z:V.z.toString()});
}

console.logS = function(V)
{
  console.log({x:V.x.toString(),y:V.y.toString(),z:V.z.toString(),vx:V.vx.toString(),vy:V.vy.toString(),vz:V.vz.toString()});
}
/*
*/

var nDigits = 22;

var view = 3;  // Initial view (z-axis)
var init = 0;  // Initial set (circular)
var corr = 1;  // Tr to do 2nd order corrections?
var doCSV = 0;

// Force modifiers:
// m - Yukawa mass (inverse range) in 1/AU
// y - Yukawa coupling constant
// r0 - Cubic galileon length scale
// a0 - MOND acceleration scale
var MOD;

{
  const params = new URLSearchParams(window.location.search);
  let m = 0;
  let y = 0;
  let r0 = 0;
  let a0 = 0;
  if (params.get('view')) view = params.get('view');
  if (params.get('init')) init = 1*params.get('init');
  if (params.get('corr')) corr = 1*params.get('corr');
  if (params.get('m')) m = 1*params.get('m'); // Yukawa mass (inverse range)
  if (params.get('y')) y = 1*params.get('y'); // Yukawa coupling constant
  if (params.get('r0')) r0 = 1*params.get('r0'); // Cubic galileon range scale
  if (params.get('a0')) a0 = 1*params.get('a0'); // MOND acceleration scale
  if (params.get('doCSV')) doCSV = params.get('doCSV');
  if (params.get('digits')) nDigits = parseInt(params.get('digits'));

  MOD = { w: new Decimal(0), m: new Decimal(m), y: new Decimal(y), r0: new Decimal(r0), a0: new Decimal(a0) };
}

Decimal.set({ precision: nDigits });
const epsilon = new Decimal("1e-" + nDigits.toString());

var strLog = "";

const GM = 1.32712440018e2; // Mm^3/kg/s^2
const GMD = new Decimal(GM);
const AU = 1.495978707e5;   // Mm

var time = 0;

var traceT = new Decimal(0);   // Inertial system
var traceW = new Decimal(0);   // Derotated satellite-fixed system

var stdevT = new Decimal(0);
var stdevW = new Decimal(0);

var Tsum = new Decimal(0);
var Tdev = new Decimal(0);
var Tnum = new Decimal(0);
var Wsum = new Decimal(0);
var Wdev = new Decimal(0);
var Wnum = new Decimal(0);


const camera = { x: 0, y: 0, z: 10, phi: 0, theta: 0 };
var spacecrafts;
var states;
var refstate;
var scale;
var step;
var indices;

class State
{
  constructor(x, y, z, vx, vy, vz, r, g, b)
  {
    this.x = new Decimal(x);
    this.y = new Decimal(y);
    this.z = new Decimal(z);
    this.vx = new Decimal(vx);
    this.vy = new Decimal(vy);
    this.vz = new Decimal(vz);
    this.r = r;
    this.g = g;
    this.b = b;
  }

  add(otherState)
  {
    let newX = this.x.plus(otherState.x);
    let newY = this.y.plus(otherState.y);
    let newZ = this.z.plus(otherState.z);
    let newVx = this.vx.plus(otherState.vx);
    let newVy = this.vy.plus(otherState.vy);
    let newVz = this.vz.plus(otherState.vz);

    return new State(newX, newY, newZ, newVx, newVy, newVz, this.r, this.g, this.b);
  }

  multiply(dt)
  {
    //return new State(this.x*dt, this.y*dt, this.z*dt, this.vx*dt, this.vy*dt, this.vz*dt, this.r, this.g, this.b);
    return new State(this.x.times(dt), this.y.times(dt), this.z.times(dt), this.vx.times(dt), this.vy.times(dt), this.vz.times(dt), this.r, this.g, this.b);
  }
}

//const nullstate = {x: 0, y: 0, z: 0};
const nullstate = new State(0, 0, 0, 0, 0, 0, 0, 0, 0);

const spacecraftsets = [
  {
    name: "Circular",
    tail: 600,
    scale: AU,
    step: 600,
    refstate: new State(AU, 0, 0, 0, 30e-3, 0),
    states: [
      new State(-1, -1, -3, 1e-8, 30e-3, 2e-7, 255, 128, 128),
      new State(+1, -1, 2, 0, 30e-3-4e-7, 0, 128, 255, 128),
      new State(-1, 1, 5, 0, 30e-3, 2e-7, 128, 128, 255),
      new State(+1, 2, 2, 0, 30e-3-4e-7, 0, 255, 255, 0)
    ]
  },
  {
    name: "Eccentric",
    tail: 600,
    scale: AU,
    step: 600,
    refstate: new State(AU, 0, 0, 0, 35e-3, 0),
    states: [
      new State(-1, -1, -1, 0e-7, 35e-3-15e-8, 1e-7, 255, 128, 128),
      new State(+1, -1, 1, 0e-7, 35e-3-50e-8, 0, 128, 255, 128),
      new State(-1, 1, -2, 0e-7, 35e-3-15e-8, 1e-7, 128, 128, 255),
      new State(+1, 1, 3, 0e-7, 35e-3-50e-8, 0e-7, 255, 255, 0)
    ]
  },
  {
    name: "High ecc.",
    tail: 600,
    scale: AU,
    step: 600,
    refstate: new State(0.6 * AU, 0, 0, 0, 0.0485, 0),
    states: [
      new State(-0.5, -0.5, -0.5, 1e-7, 0.0485+1.7e-7, 5e-7, 255, 128, 128),
      new State(+0.5, -0.5, 0.5, 1e-7, 0.0485-1.7e-7, -8e-7, 128, 255, 128),
      new State(-0.5, 0.5, 0.5, 1e-7, 0.0485+1.7e-7, -4e-7, 128, 128, 255),
      new State(+0.5, 0.5, -0.5, -1e-7, 0.0485-1.7e-7, 3e-7, 255, 255, 0)
    ]
  },
  {
    name: "High ecc. (tiny)",
    tail: 600,
    scale: AU,
    step: 30,
    refstate: new State(0.6 * AU, 0, 0, 0, 0.0485, 0),
    states: [
      new State(-0.5e-1, -0.5e-1, -0.5e-1, 1e-8, 0.0485+1.7e-8, 5e-8, 255, 128, 128),
      new State(+0.5e-1, -0.5e-1, 0.5e-1, 1e-8, 0.0485-1.7e-8, -8e-8, 128, 255, 128),
      new State(-0.5e-1, 0.5e-1, 0.5e-1, 1e-8, 0.0485+1.7e-8, -4e-8, 128, 128, 255),
      new State(+0.5e-1, 0.5e-1, -0.5e-1, -1e-8, 0.0485-1.7e-8, 3e-8, 255, 255, 0)
    ]
  },
  {
    name: "High ecc. (large)",
    tail: 600,
    scale: AU,
    step: 30,
    refstate: new State(0.6 * AU, 0, 0, 0, 0.0485, 0),
    states: [
      new State(-0.5e1, -0.5e1, -0.5e1, 2e-6, 0.0485+1.8e-6, -2.2e-6, 255, 128, 128),
      new State(+0.5e1, -0.5e1, 0.5e1, 1e-6, 0.0485-0.7e-6, 2.2e-6, 128, 255, 128),
      new State(-0.5e1, 0.5e1, 0.5e1, -1e-6, 0.0485+0.7e-6, 2.2e-6, 128, 128, 255),
      new State(+0.5e1, 0.5e1, -0.5e1, -2e-6, 0.0485-1.8e-6, -2.2e-6, 255, 255, 0)
    ]
  },
  {
    name: "Medium eccentric",
    tail: 600,
    scale: 2*AU,
    step: 3600,
    refstate: new State(1*AU, 0, 0, 0, 35e-3, 0),
    states: [
      new State(-1, -1, -3, 1e-10, 35e-3, 1e-9, 255, 128, 128),
      new State(+1, -1, 2, 0, 35e-3-2e-9, 0, 128, 255, 128),
      new State(-1, 1, 5, 0, 35e-3, 1e-9, 128, 128, 255),
      new State(+1, 2, 2, 0, 35e-3-2e-9, 0, 255, 255, 0)
    ]
  },
  {
    name: "Big eccentric",
    tail: 3000,
    scale: 10*AU,
    step: 21600,
    refstate: new State(10*AU, 0, 0, 0, 8e-3, 0),
    states: [
      new State(-1, -1, -3, 1e-9, 8e-3,      1e-8, 255, 128, 128),
      new State(+1, -1,  2, 0,    8e-3-2e-8, 0,    128, 255, 128),
      new State(-1,  1,  5, 0,    8e-3,      1e-8, 128, 128, 255),
      new State(+1,  2,  2, 0,    8e-3-2e-8, 0,    255, 255, 0)
    ]
  },
  {
    name: "Gigantic eccentric",
    tail: 3000,
    scale: 30*AU,
    step: 2*86400,
    refstate: new State(30*AU, 0, 0, 0, 5e-3, 0),
    states: [
      new State(-1, -1, -3, 1e-9, 5e-3,      1e-10, 255, 128, 128),
      new State(+1, -1,  2, 0,    5e-3-2e-9, 0,     128, 255, 128),
      new State(-1,  1,  5, 0,    5e-3,      1e-10, 128, 128, 255),
      new State(+1,  2,  2, 0,    5e-3-2e-9, 0,     255, 255, 0)
    ]
  },
  {
    name: "Tiny eccentric",
    tail: 300,
    scale: 0.1*AU,
    step: 90,
    refstate: new State(0.1*AU, 0, 0, 0, 100e-3, 0),
    states: [
      new State(-0.1, -0.1, -0.2,  2e-8, 100e-3,      0e-8, 255, 128, 128),
      new State(+0.1, -0.1,  0.1, -1e-8, 100e-3-0e-8, 0,    128, 255, 128),
      new State(-0.1,  0.1,  0.2,  1e-8, 100e-3,      0e-8, 128, 128, 255),
      new State(+0.1, -0.2,  0.1, -3e-8, 100e-3+0e-5, 0,    255, 255, 0)
    ]
  },
  {
    name: "test",
    tail: 300,
    scale: AU,
    step: 6,
    refstate: new State(AU, 0, 0, 0, 30e-3, 0),
    states: [
      new State( 0, 0, 0, 0, 30e-3, 0,    255, 128, 128),
      new State( 0, 1, 0, 0, 30e-3, 0,    128, 255, 128),
      new State( 0, 0, 1, 0, 30e-3, 0,    128, 128, 255),
      new State( 1, 0, 0, 0, 30e-3, 0,    255, 255, 0)
    ]
  },
  {
    name: "High ecc. [alt]",
    tail: 600,
    scale: AU,
    step: 600,
    refstate: new State(0.6 * AU, 0, 0, 0, 0.0485, 0),
    states: [
      new State(-0.5, -0.5, -0.5,  1e-7, 0.0485+1.6e-7,  5e-7, 255, 128, 128),
      new State(+0.5, -0.5,  0.5,  1e-7, 0.0485-1.8e-7, -8e-7, 128, 255, 128),
      new State(-0.5,  0.5,  0.5,  1e-7, 0.0485+1.6e-7, -4e-7, 128, 128, 255),
      new State(+0.5,  0.5, -0.5, -1e-7, 0.0485-1.8e-7,  3e-7, 255, 255, 0)
    ]
  },
];

function saveAll()
{
  let savedata = JSON.stringify(
  {
    init: init,
    view: view,
    corr: corr,
    sgna: 1,
    lina: 0,
    states: states,
    refstate: refstate,
    scale: scale,
    step: step,
    time: time,
    MOD: MOD,
    traceT: traceT,
    traceW: traceW,
    stdevT: stdevT,
    stdevW: stdevW,
    Tsum: Tsum,
    Tdev: Tdev,
    Tnum: Tnum,
    Wsum: Wsum,
    Wdev: Wdev,
    Wnum: Wnum,
    top: 
    {
      size: topView.bufferSize,
      buffer: topView.buffer
    },
    rts: {M: rtsM, T: rtsT, W: rtsW},
    camera: camera
  });

  let blob = new Blob([savedata], {type: "application/json"});
  let url = URL.createObjectURL(blob);
  let a = document.createElement("a");
  a.href = url;
  a.download=document.title + ".json";
  a.click();
}

function loadAll()
{
  stop();

  let fileInput = document.createElement("input");
  fileInput.type = "file";
  fileInput.accept = ".json, application/json";
  fileInput.addEventListener("change", () =>
  {
    let file = fileInput.files[0];
    fileName = file.name.replace(/.json$/,"");
    let reader = new FileReader();
    reader.addEventListener("load", () =>
    {
      document.getElementById("init").disabled = true;

      let data = JSON.parse(reader.result);
      init = data.init;
      view = data.view;
      corr = data.corr;
      scale = data.scale;
      step = data.step;
      time = parseFloat(data.time);
      MOD = data.MOD;
       MOD.w = new Decimal(MOD.w);
       MOD.m = new Decimal(MOD.m);
       MOD.y = new Decimal(MOD.y);
      traceT = new Decimal(data.traceT);
      traceW = new Decimal(data.traceW);
      stdevT = new Decimal(data.stdevT);
      stdevW = new Decimal(data.stdevW);
      Tsum = new Decimal(data.Tsum);
      Tdev = new Decimal(data.Tdev);
      Tnum = data.Tnum;
      Wsum = new Decimal(data.Wsum);
      Wdev = new Decimal(data.Wdev);
      Wnum = data.Wnum;
      topView.bufferSize = data.top.size;
      topView.buffer = data.top.buffer;
      rtsM = [];
      rtsT = [];
      rtsW = [];
      for (let K = 0; K < 4; K++)
      {
        function doData(d)
        {
          const dD = [];
          for (let i = 0; i < d.length; i++)
          {
            if ("ex" in d[i])
            {
              dD.push({
               r: new Decimal(d[i].r),
               n: {x: new Decimal(d[i].n.x), y: new Decimal(d[i].n.y), z: new Decimal(d[i].n.z)},
               r12: {x: new Decimal(d[i].r12.x), y: new Decimal(d[i].r12.y), z: new Decimal(d[i].r12.z)},
               r13: {x: new Decimal(d[i].r13.x), y: new Decimal(d[i].r13.y), z: new Decimal(d[i].r13.z)},
               r14: {x: new Decimal(d[i].r14.x), y: new Decimal(d[i].r14.y), z: new Decimal(d[i].r14.z)},
               ex: {x: new Decimal(d[i].ex.x), y: new Decimal(d[i].ex.y), z: new Decimal(d[i].ex.z)},
               ey: {x: new Decimal(d[i].ey.x), y: new Decimal(d[i].ey.y), z: new Decimal(d[i].ey.z)},
               ez: {x: new Decimal(d[i].ez.x), y: new Decimal(d[i].ez.y), z: new Decimal(d[i].ez.z)}});
            }
            else if ("r1" in d[i])
            {
              dD.push({
               r: new Decimal(d[i].r),
               n: {x: new Decimal(d[i].n.x), y: new Decimal(d[i].n.y), z: new Decimal(d[i].n.z)},
               r12: {x: new Decimal(d[i].r12.x), y: new Decimal(d[i].r12.y), z: new Decimal(d[i].r12.z)},
               r13: {x: new Decimal(d[i].r13.x), y: new Decimal(d[i].r13.y), z: new Decimal(d[i].r13.z)},
               r14: {x: new Decimal(d[i].r14.x), y: new Decimal(d[i].r14.y), z: new Decimal(d[i].r14.z)},
               r1: {x: new Decimal(d[i].r1.x), y: new Decimal(d[i].r1.y), z: new Decimal(d[i].r1.z)}});
            }
            else
            {
              dD.push({
               r: new Decimal(d[i].r),
               n: {x: new Decimal(d[i].n.x), y: new Decimal(d[i].n.y), z: new Decimal(d[i].n.z)},
               r12: {x: new Decimal(d[i].r12.x), y: new Decimal(d[i].r12.y), z: new Decimal(d[i].r12.z)},
               r13: {x: new Decimal(d[i].r13.x), y: new Decimal(d[i].r13.y), z: new Decimal(d[i].r13.z)},
               r14: {x: new Decimal(d[i].r14.x), y: new Decimal(d[i].r14.y), z: new Decimal(d[i].r14.z)}});
            }
          }
          return dD;
        }
        rtsM.push(new TS(data.rts.M[K].size));
        rtsM[K].data = doData(data.rts.M[K].data);
        rtsT.push(new TS(data.rts.T[K].size));
        rtsT[K].data = doData(data.rts.T[K].data);
        rtsW.push(new TS(data.rts.W[K].size));
        rtsW[K].data = doData(data.rts.W[K].data);
      }
      camera.x = data.camera.x;
      camera.y = data.camera.y;
      camera.z = data.camera.z;

      document.getElementById("init").value = init;
      doResize();
      onInit(init);

      // Must be done after onInit();
      refstate = new State(data.refstate.x, data.refstate.y, data.refstate.z,
                           data.refstate.vx, data.refstate.vy, data.refstate.vz,
                           data.refstate.r, data.refstate.g, data.refstate.b);
      states = [];
      data.states.forEach(s =>
      {
        states.push(new State(s.x, s.y, s.z, s.vx, s.vy, s.vz, s.r, s.g, s.b));
      });
      document.getElementById("view").value = view;
      setView(view);
      fileInput.remove();
    });
    reader.readAsText(file);
  });
  fileInput.click();
}

function stormerRichardson(s, dt, refstate)
{
  function a(s, refstate)
  {
    var r = Decimal.sqrt((s.x.plus(refstate.x)).pow(2).plus((s.y.plus(refstate.y)).pow(2)).plus((s.z.plus(refstate.z)).pow(2)));
    var r3 = r.pow(3);
    var GMY = GMD.times(new Decimal(1).plus(MOD.y.times(new Decimal(1).minus(new Decimal(1).plus(r.times(MOD.m))).times(r.neg().times(MOD.m).exp()))));
    if (!MOD.r0.isZero()) GMY = GMY.times(new Decimal(1).plus(Decimal.pow(r.div(MOD.r0), new Decimal(1.5))));
    if (!MOD.a0.isZero()) GMY = GMY.times(Decimal.sqrt(new Decimal(0.5).plus(Decimal.sqrt(new Decimal(0.25).plus(MOD.a0.times(r).times(r).div(GMD).pow(2))))));


//    if (!MOD.a0.isZero()) GMY = GMY.times(Decimal.sqrt(new Decimal(1).div(new Decimal(1).plus(MOD.a0.times(MOD.a0).times(r3).times(r).div(GMY).div(GMY)))));

    return new State(new Decimal(0), new Decimal(0), new Decimal(0),
                     GMY.neg().times(s.x.plus(refstate.x)).div(r3),
                     GMY.neg().times(s.y.plus(refstate.y)).div(r3),
                     GMY.neg().times(s.z.plus(refstate.z)).div(r3));
  }

  function leapfrogStep(s, dt, refstate)
  {
    // Half-step velocity update
    var a0 = a(s, refstate);
    var vx_half = s.vx.plus(dt.div(2).times(a0.vx));
    var vy_half = s.vy.plus(dt.div(2).times(a0.vy));
    var vz_half = s.vz.plus(dt.div(2).times(a0.vz));

    // Full-step position update
    var x_new = s.x.plus(dt.times(vx_half));
    var y_new = s.y.plus(dt.times(vy_half));
    var z_new = s.z.plus(dt.times(vz_half));

    // Full-step acceleration update
    var s_new = new State(x_new, y_new, z_new, vx_half, vy_half, vz_half, s.r, s.g, s.b);
    var a_new = a(s_new, refstate);

    // Half-step velocity update
    var vx_new = vx_half.plus(dt.div(2).times(a_new.vx));
    var vy_new = vy_half.plus(dt.div(2).times(a_new.vy));
    var vz_new = vz_half.plus(dt.div(2).times(a_new.vz));

    return new State(x_new, y_new, z_new, vx_new, vy_new, vz_new, s.r, s.g, s.b);
  }

  // Perform two leapfrog steps with half the time step
  var halfDt = dt.div(2);
  var s_half1 = leapfrogStep(s, halfDt, refstate);
  var s_half2 = leapfrogStep(s_half1, halfDt, refstate);

  // Perform one leapfrog step with the full time step
  var s_full = leapfrogStep(s, dt, refstate);

  // Richardson extrapolation
  var x_new = s_half2.x.times(4).minus(s_full.x).div(3);
  var y_new = s_half2.y.times(4).minus(s_full.y).div(3);
  var z_new = s_half2.z.times(4).minus(s_full.z).div(3);
  var vx_new = s_half2.vx.times(4).minus(s_full.vx).div(3);
  var vy_new = s_half2.vy.times(4).minus(s_full.vy).div(3);
  var vz_new = s_half2.vz.times(4).minus(s_full.vz).div(3);

  return new State(x_new, y_new, z_new, vx_new, vy_new, vz_new, s.r, s.g, s.b);
}

const integrator = stormerRichardson;

var sunZ;

function transformCoordinates(stateVectors, refstate)
{
  // Step 1: Find the geometric center
  let centerX = new Decimal(0), centerY = new Decimal(0), centerZ = new Decimal(0);
  for (let i = 0; i < stateVectors.length; i++)
  {
    centerX = centerX.plus(stateVectors[i].x);
    centerY = centerY.plus(stateVectors[i].y);
    centerZ = centerZ.plus(stateVectors[i].z);
  }
  centerX = centerX.div(stateVectors.length);
  centerY = centerY.div(stateVectors.length);
  centerZ = centerZ.div(stateVectors.length);

  // Step 2: Designate the new origin O'
  let newOrigin = {x: centerX, y: centerY, z: centerZ};

  // Step 3: Establish the new coordinate system
  let zPrimeAxis = {x: centerX.plus(refstate.x), y: centerY.plus(refstate.y), z: centerZ.plus(refstate.z)};
  let zPrimeAxisLength = Decimal.sqrt(zPrimeAxis.x.pow(2).plus(zPrimeAxis.y.pow(2)).plus(zPrimeAxis.z.pow(2)));
  zPrimeAxis.x = zPrimeAxis.x.div(zPrimeAxisLength);
  zPrimeAxis.y = zPrimeAxis.y.div(zPrimeAxisLength);
  zPrimeAxis.z = zPrimeAxis.z.div(zPrimeAxisLength);

  let xPrimeAxis = {x: zPrimeAxis.y.neg(), y: zPrimeAxis.x, z: new Decimal(0)};
  let xPrimeAxisLength = Decimal.sqrt(xPrimeAxis.x.pow(2).plus(xPrimeAxis.y.pow(2)).plus(xPrimeAxis.z.pow(2)));
  xPrimeAxis.x = xPrimeAxis.x.div(xPrimeAxisLength);
  xPrimeAxis.y = xPrimeAxis.y.div(xPrimeAxisLength);
  xPrimeAxis.z = xPrimeAxis.z.div(xPrimeAxisLength);

  // Step 4: Express the coordinates of each spacecraft in the new coordinate system
  let transformedCoordinates = [];
  for (let i = 0; i < stateVectors.length; i++) {
    let dx = stateVectors[i].x.minus(newOrigin.x);
    let dy = stateVectors[i].y.minus(newOrigin.y);
    let dz = stateVectors[i].z.minus(newOrigin.z);

    let xPrime = dx.times(xPrimeAxis.x).plus(dy.times(xPrimeAxis.y)).plus(dz.times(xPrimeAxis.z));
    let yPrime = dx.times(zPrimeAxis.y.times(xPrimeAxis.z).minus(zPrimeAxis.z.times(xPrimeAxis.y)))
               .plus(dy.times(zPrimeAxis.z.times(xPrimeAxis.x).minus(zPrimeAxis.x.times(xPrimeAxis.z))))
               .plus(dz.times(zPrimeAxis.x.times(xPrimeAxis.y).minus(zPrimeAxis.y.times(xPrimeAxis.x))));
    let zPrime = dx.times(zPrimeAxis.x).plus(dy.times(zPrimeAxis.y)).plus(dz.times(zPrimeAxis.z));

    transformedCoordinates.push({x: xPrime, y: yPrime, z: zPrime,
            vx: stateVectors[i].vx, vy: stateVectors[i].vy, vz: stateVectors[i].vz,
            r: stateVectors[i].r, g: stateVectors[i].g, b: stateVectors[i].b});
  }

  sunZ = norm(refstate).neg();  // Approximate position of the Sun

  return transformedCoordinates;
}

function determinant(r1, r2, r3, r4)
{
  return r2.x.times(r3.y).times(r4.z).minus(
         r1.x.times(r3.y).times(r4.z)).minus(
         r2.y.times(r3.x).times(r4.z)).plus(
         r1.y.times(r3.x).times(r4.z)).plus(
         r1.x.times(r2.y).times(r4.z)).minus(
         r1.y.times(r2.x).times(r4.z)).minus(
         r2.x.times(r3.z).times(r4.y)).plus(
         r1.x.times(r3.z).times(r4.y)).plus(
         r2.z.times(r3.x).times(r4.y)).minus(
         r1.z.times(r3.x).times(r4.y)).minus(
         r1.x.times(r2.z).times(r4.y)).plus(
         r1.z.times(r2.x).times(r4.y)).plus(
         r2.y.times(r3.z).times(r4.x)).minus(
         r1.y.times(r3.z).times(r4.x)).minus(
         r2.z.times(r3.y).times(r4.x)).plus(
         r1.z.times(r3.y).times(r4.x)).plus(
         r1.y.times(r2.z).times(r4.x)).minus(
         r1.z.times(r2.y).times(r4.x)).minus(
         r1.x.times(r2.y).times(r3.z)).plus(
         r1.y.times(r2.x).times(r3.z)).plus(
         r1.x.times(r2.z).times(r3.y)).minus(
         r1.z.times(r2.x).times(r3.y)).minus(
         r1.y.times(r2.z).times(r3.x)).plus(
         r1.z.times(r2.y).times(r3.x));
}

function volume(spacecrafts)
{
  return Decimal.div(determinant(spacecrafts[0], spacecrafts[1], spacecrafts[2], spacecrafts[3]), 6);
}

function area(s1, s2, s3)
{
  var d = Decimal.pow(s2.z.minus(s3.z), 2).plus(Decimal.pow(s2.y.minus(s3.y), 2)).plus(Decimal.pow(s2.x.minus(s3.x), 2))
          .times(Decimal.pow(s1.x.minus(s3.x), 2).plus(Decimal.pow(s1.z.minus(s3.z), 2)).plus(Decimal.pow(s1.y.minus(s3.y), 2)));
  d = d.minus(Decimal.pow(s2.x.minus(s3.x).times(s1.x.minus(s3.x)).plus(s1.z.minus(s3.z).times(s2.z.minus(s3.z))).plus(s1.y.minus(s3.y).times(s2.y.minus(s3.y))), 2));
  return d.sqrt();
}

function L(r1,r2)
{
  return Decimal.pow(r1.x.minus(r2.x), 2).plus(Decimal.pow(r1.y.minus(r2.y), 2)).plus(Decimal.pow(r1.z.minus(r2.z), 2)).sqrt();
}

function X(u,v)
{
  return {x:u.y.times(v.z).minus(u.z.times(v.y)), y:u.z.times(v.x).minus(u.x.times(v.z)), z:u.x.times(v.y).minus(u.y.times(v.x))};
}

function dot(u,v)
{
  return u.x.times(v.x).plus(u.y.times(v.y)).plus(u.z.times(v.z));
}

function norm(R)
{
  return Decimal.pow(R.x, 2).plus(Decimal.pow(R.y, 2)).plus(Decimal.pow(R.z, 2)).sqrt();
}

function smul(s, v)
{
  return {x:s.times(v.x), y:s.times(v.y), z:s.times(v.z)};
}

function Mmul(M, v)
{
  return {x:M[0][0].times(v.x).plus(M[0][1].times(v.y)).plus(M[0][2].times(v.z)),
          y:M[1][0].times(v.x).plus(M[1][1].times(v.y)).plus(M[1][2].times(v.z)),
          z:M[2][0].times(v.x).plus(M[2][1].times(v.y)).plus(M[2][2].times(v.z))};
}

function MP(M, K)
{
  let R = [[new Decimal(0), new Decimal(0), new Decimal(0)], 
           [new Decimal(0), new Decimal(0), new Decimal(0)], 
           [new Decimal(0), new Decimal(0), new Decimal(0)]];

  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
  {
    R[i][j] = M[i][j].plus(K[i][j]);
  }
  return R;
}

function MM(M, K)
{
  let R = [[new Decimal(0), new Decimal(0), new Decimal(0)], 
           [new Decimal(0), new Decimal(0), new Decimal(0)], 
           [new Decimal(0), new Decimal(0), new Decimal(0)]];

  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
      for (let k = 0; k < 3; k++)
  {
    R[i][j] = R[i][j].plus(M[i][k].times(K[k][j]));
  }
  return R;
}

function sM(s, M)
{
  let R = [[new Decimal(0), new Decimal(0), new Decimal(0)], 
           [new Decimal(0), new Decimal(0), new Decimal(0)], 
           [new Decimal(0), new Decimal(0), new Decimal(0)]];

  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
  {
    R[i][j] = s.times(M[i][j]);
  }
  return R;
}

function vadd(u, v)
{
  return {x:u.x.plus(v.x), y:u.y.plus(v.y), z:u.z.plus(v.z)};
}

function inv3x3M(M)
{
  const a = M[0][0], b = M[0][1], c = M[0][2];
  const d = M[1][0], e = M[1][1], f = M[1][2];
  const g = M[2][0], h = M[2][1], i = M[2][2];

  const D = a.times(e.times(i).minus(f.times(h))).minus(b.times(d.times(i).minus(f.times(g)))).plus(c.times(d.times(h).minus(e.times(g))));

  return [[e.times(i).minus(f.times(h)).div(D), c.times(h).minus(b.times(i)).div(D), b.times(f).minus(c.times(e)).div(D)],
          [f.times(g).minus(d.times(i)).div(D), a.times(i).minus(c.times(g)).div(D), c.times(d).minus(a.times(f)).div(D)],
          [d.times(h).minus(e.times(g)).div(D), b.times(g).minus(a.times(h)).div(D), a.times(e).minus(b.times(d)).div(D)]];
}

class TS
{
  constructor(size)
  {
    this.size = size;
    this.data = [];
  }
  push(element)
  {
    this.data.push(element);
    if (this.data.length > this.size) this.data.shift();
  }
}

var rtsM = [];
var rtsT = [];
var rtsW = [];

for (let K = 0; K < 4; K++)
{
  rtsM.push(new TS(5));
  rtsT.push(new TS(5));
  rtsW.push(new TS(5));
}

function trT(states, DT, K=0)
{
  let r1 = {x:states[K].x.plus(refstate.x), y:states[K].y.plus(refstate.y), z: states[K].z.plus(refstate.z)};
  let r12 = {x:states[(K+1)%4].x.minus(states[K].x), y:states[(K+1)%4].y.minus(states[K].y), z: states[(K+1)%4].z.minus(states[K].z)};
  let r13 = {x:states[(K+2)%4].x.minus(states[K].x), y:states[(K+2)%4].y.minus(states[K].y), z: states[(K+2)%4].z.minus(states[K].z)};
  let r14 = {x:states[(K+3)%4].x.minus(states[K].x), y:states[(K+3)%4].y.minus(states[K].y), z: states[(K+3)%4].z.minus(states[K].z)};

  let n = refstate;
  let r = norm(n);
  n = smul(r.pow(-1), n);

  rtsT[K].push({r1: r1, r12:r12, r13:r13, r14:r14, r:r, n:n});

  return doTR(rtsT[K], DT);
}

function trM(states, DT, K=0)
{
  // First we obtain the six edges of the tetrahedron
  let a = L(states[K], states[(K+1)%4]);
  let b = L(states[K], states[(K+2)%4]);
  let c = L(states[(K+1)%4], states[(K+2)%4]);
  let d = L(states[K], states[(K+3)%4]);
  let e = L(states[(K+1)%4], states[(K+3)%4]);
  let f = L(states[(K+2)%4], states[(K+3)%4]);

  // Next, we establish the coordinates of the four vertices
  // in a reference frame affixed to the satellite constellation
  let r12={x:a, y:new Decimal(0), z:new Decimal(0)};
  let r13={x:a.pow(2).plus(b.pow(2)).minus(c.pow(2)).div(a.times(2)), 
           y:Decimal.sqrt(b.pow(2).minus(a.pow(2).plus(b.pow(2)).minus(c.pow(2)).div(a.times(2)).pow(2))), 
           z:new Decimal(0)};
  let r14={x:a.pow(2).plus(d.pow(2)).minus(e.pow(2)).div(a.times(2)),
           y:b.pow(2).plus(d.pow(2)).minus(f.pow(2)).minus(a.pow(2).plus(d.pow(2)).minus(e.pow(2)).times(a.pow(2).plus(b.pow(2)).minus(c.pow(2))).div(a.pow(2).times(2)))
             .div(Decimal.sqrt(b.pow(2).minus(a.pow(2).plus(b.pow(2)).minus(c.pow(2)).div(a.times(2)).pow(2))).times(2)),
           z:Decimal.sqrt(d.pow(2).minus(a.pow(2).plus(d.pow(2)).minus(e.pow(2)).div(a.times(2)).pow(2))
             .minus(b.pow(2).plus(d.pow(2)).minus(f.pow(2)).minus(a.pow(2).plus(d.pow(2)).minus(e.pow(2)).times(a.pow(2).plus(b.pow(2)).minus(c.pow(2))).div(a.pow(2).times(2)))
             .div(Decimal.sqrt(b.pow(2).minus(a.pow(2).plus(b.pow(2)).minus(c.pow(2)).div(a.times(2)).pow(2))).times(2)).pow(2)))};

  // We now compute the unit vectors of the satellite-fixed frame
  let ex = vadd(states[(K+1)%4], smul(new Decimal(-1), states[K]));
  let ey = vadd(states[(K+2)%4], smul(new Decimal(-1), states[K]));
  let ez = X(ex, ey);
  ex = smul(norm(ex).pow(-1), ex);
  ez = smul(norm(ez).pow(-1), ez);
  ey = X(ez, ex);

  // The Sun-to-origin direction in this frame is obtained by
  // projecting that vector onto the unit vectors
  let n = refstate;
  let r = norm(n);
  n = {x: dot(n, ex), y: dot(n, ey), z: dot(n, ez)};
  n = smul(norm(n).pow(-1), n);

  // There is a "chiral ambiguity" that we resolve by presuming
  // that the constellation is smart enough to distinguish itself
  // from its mirror image. In our case, we check if the handedness
  // of the tetrahedron matches that of the reference frame, by
  // projecting vertex D onto the ez axis:
  if (dot(ez,vadd(states[(K+3)%4],smul(new Decimal(-1),states[K]))).times(r14.z).lt(0)) r14.z = r14.z.neg();

  rtsM[K].push({r12:r12, r13:r13, r14:r14, r:r, n:n, ex:ex, ey:ey, ez:ez});

  return new Decimal(0);
}

function SagnacDiff(ri, rj, vi, vj)
{
  const c = new Decimal("299.792458"); // Mm/s
  const c2 = c.pow(2);

  let t1i = dot(vi,ri).plus(Decimal.sqrt(dot(vi,ri).pow(2).plus(c2.minus(dot(vi,vi)).times(dot(ri,ri))))).div(c2.minus(dot(vi,vi)));
  let rij = vadd(rj, smul(new Decimal(-1), ri));
  let vij = vadd(vj, smul(new Decimal(-1), vi));
  let bij = dot(vadd(rij, smul(t1i, vij)), vj);
  let tij = bij.plus(Decimal.sqrt(bij.pow(2).plus(c2.minus(dot(vj,vj)).times(dot(vadd(rij, smul(t1i, vij)),vadd(rij, smul(t1i, vij))))))).div(c2.minus(dot(vj,vj)));
  let tj1 = Decimal.sqrt(dot(vadd(rj, smul(t1i.plus(tij), vj)), vadd(rj, smul(t1i.plus(tij), vj)))).div(c);

  let t1j = dot(vj,rj).plus(Decimal.sqrt(dot(vj,rj).pow(2).plus(c2.minus(dot(vj,vj)).times(dot(rj,rj))))).div(c2.minus(dot(vj,vj)));
  let rji = vadd(ri, smul(new Decimal(-1), rj));
  let vji = vadd(vi, smul(new Decimal(-1), vj));
  let bji = dot(vadd(rji, smul(t1j, vji)), vi);
  let tji = bji.plus(Decimal.sqrt(bji.pow(2).plus(c2.minus(dot(vi,vi)).times(dot(vadd(rji, smul(t1j, vji)),vadd(rji, smul(t1j, vji))))))).div(c2.minus(dot(vi,vi)));
  let ti1 = Decimal.sqrt(dot(vadd(ri, smul(t1j.plus(tji), vi)), vadd(ri, smul(t1j.plus(tji), vi)))).div(c);

  let dt = t1i.minus(ti1);
  dt = dt.plus(tij.minus(tji));
  dt = dt.plus(tj1.minus(t1j));

  return dt;
}

function trW(DT, K=0)
{
  // trW is computed by correcting trM for rotation.
  // The rotation is estimated from a "measurement" of the Sagnac-type
  // observable. To simulate this measurement, we first need to
  // compute that observable. We do so for the middle of a triplet of
  // position observables, as this will be the midpoint used later
  // on for numerically calculating acceleration.
  if (rtsT[K].data.length != 5) return new Decimal(NaN);

  // Midpoint positions
  let r12 = rtsT[K].data[2].r12;
  let r13 = rtsT[K].data[2].r13;
  let r14 = rtsT[K].data[2].r14;

  // Velocities of satellites 234 wrt. satellite 1

  let v12 = smul(new Decimal(1).div(DT.times(12)), vadd(vadd(rtsT[K].data[0].r12, smul(new Decimal(-8), rtsT[K].data[1].r12)), vadd(smul(new Decimal(8), rtsT[K].data[3].r12), smul(new Decimal(-1), rtsT[K].data[4].r12))));
  let v13 = smul(new Decimal(1).div(DT.times(12)), vadd(vadd(rtsT[K].data[0].r13, smul(new Decimal(-8), rtsT[K].data[1].r13)), vadd(smul(new Decimal(8), rtsT[K].data[3].r13), smul(new Decimal(-1), rtsT[K].data[4].r13))));
  let v14 = smul(new Decimal(1).div(DT.times(12)), vadd(vadd(rtsT[K].data[0].r14, smul(new Decimal(-8), rtsT[K].data[1].r14)), vadd(smul(new Decimal(8), rtsT[K].data[3].r14), smul(new Decimal(-1), rtsT[K].data[4].r14))));

  // The actual Sagnac-type time differences
  let w123 = SagnacDiff(r12, r13, v12, v13);
  let w134 = SagnacDiff(r13, r14, v13, v14);
  let w142 = SagnacDiff(r14, r12, v14, v12);

  // From this point on, all we are allowed to use are w1..w3 and rtsM
  // as only this information will be available on board. We are now
  // in the satellite-fixed noninertial reference frame.
  //

  r12 = rtsM[K].data[2].r12;
  r13 = rtsM[K].data[2].r13;
  r14 = rtsM[K].data[2].r14;

  // Velocities of satellites 234 wrt. satellite 1

v12 = smul(new Decimal(1).div(DT.times(12)), vadd(vadd(rtsM[K].data[0].r12, smul(new Decimal(-8), rtsM[K].data[1].r12)), vadd(smul(new Decimal(8), rtsM[K].data[3].r12), smul(new Decimal(-1), rtsM[K].data[4].r12))));
v13 = smul(new Decimal(1).div(DT.times(12)), vadd(vadd(rtsM[K].data[0].r13, smul(new Decimal(-8), rtsM[K].data[1].r13)), vadd(smul(new Decimal(8), rtsM[K].data[3].r13), smul(new Decimal(-1), rtsM[K].data[4].r13))));
v14 = smul(new Decimal(1).div(DT.times(12)), vadd(vadd(rtsM[K].data[0].r14, smul(new Decimal(-8), rtsM[K].data[1].r14)), vadd(smul(new Decimal(8), rtsM[K].data[3].r14), smul(new Decimal(-1), rtsM[K].data[4].r14))));

  // We now solve for the angular velocity vector w.
  let w = {x:new Decimal(0), y:new Decimal(0), z:new Decimal(0)};
  let eps = epsilon.sqrt();
  let dx = {x:eps, y:new Decimal(0), z:new Decimal(0)};
  let dy = {x:new Decimal(0), y:eps, z:new Decimal(0)};
  let dz = {x:new Decimal(0), y:new Decimal(0), z:eps};

  // Helper function to calculate modeled Sagnac-type observables
  function m1234(w, r12, r13, r14, v12, v13, v14)
  {
    let w12 = vadd(v12, X(r12, w));
    let w13 = vadd(v13, X(r13, w));
    let w14 = vadd(v14, X(r14, w));

    let m =
    {
      m123: SagnacDiff(r12, r13, w12, w13),
      m134: SagnacDiff(r13, r14, w13, w14),
      m142: SagnacDiff(r14, r12, w14, w12)
    };
    return m;
  };

  // The solution converges rapidly so 10 iterations are sufficient.
  let count = 10;
  let m;
  let T;
  let Dw = {x:new Decimal(0), y:new Decimal(0), z:new Decimal(0)};
  let Fw = {x:new Decimal(0), y:new Decimal(0), z:new Decimal(0)};
  let Da = {x:new Decimal(0), y:new Decimal(0), z:new Decimal(0)};
  let Fa = {x:new Decimal(0), y:new Decimal(0), z:new Decimal(0)};

  do
  {
    w = vadd(w, Dw);

    m = m1234(w, r12, r13, r14, v12, v13, v14);

    // Calculating a gradient matrix J in our solution space
    let mx = m1234(vadd(w, dx), r12, r13, r14, v12, v13, v14);
    let my = m1234(vadd(w, dy), r12, r13, r14, v12, v13, v14);
    let mz = m1234(vadd(w, dz), r12, r13, r14, v12, v13, v14);

    let J = [[(mx.m123.minus(m.m123)).div(eps), (my.m123.minus(m.m123)).div(eps), (mz.m123.minus(m.m123)).div(eps)],
             [(mx.m134.minus(m.m134)).div(eps), (my.m134.minus(m.m134)).div(eps), (mz.m134.minus(m.m134)).div(eps)],
             [(mx.m142.minus(m.m142)).div(eps), (my.m142.minus(m.m142)).div(eps), (mz.m142.minus(m.m142)).div(eps)]];

    Fw = {x:w123.minus(m.m123), y:w134.minus(m.m134), z:w142.minus(m.m142)};
    Dw = Mmul(inv3x3M(J), Fw);

  } while (--count > 0 && Decimal.sqrt(dot(Fw,Fw)).gt(epsilon));

  if (count < 1)
  {
    console.log("Insufficient accuracy computing w: " + Decimal.sqrt(dot(Fw,Fw)).toString());
  }

  // Constructing the rotation matrix using Rodrigues' formula
  let W = Decimal.sqrt(dot(w,w));
  if (W.eq(0)) W = new Decimal(1);
  let I = [[new Decimal(1), new Decimal(0), new Decimal(0)], [new Decimal(0), new Decimal(1), new Decimal(0)], [new Decimal(0), new Decimal(0), new Decimal(1)]];
  let Q = [[new Decimal(0), w.z.neg().div(W), w.y.div(W)], [w.z.div(W), new Decimal(0), w.x.neg().div(W)], [w.y.neg().div(W), w.x.div(W), new Decimal(0)]];

  let theta = W.times(DT);
  let R = MP(I, MP(sM(theta.sin(), Q), sM(new Decimal(1).minus(theta.cos()), MM(Q,Q))));

  // We now adjust the first two and last two position vectors to remove the
  // pseudoforce due to frame rotation.

  r12 = Mmul(R, Mmul(R, rtsM[K].data[0].r12));
  r13 = Mmul(R, Mmul(R, rtsM[K].data[0].r13));
  r14 = Mmul(R, Mmul(R, rtsM[K].data[0].r14));

  rtsW[K].push({r12:r12, r13:r13, r14:r14, r:rtsM[K].data[2].r, n:rtsM[K].data[2].n});

  r12 = Mmul(R, rtsM[K].data[1].r12);
  r13 = Mmul(R, rtsM[K].data[1].r13);
  r14 = Mmul(R, rtsM[K].data[1].r14);

  rtsW[K].push({r12:r12, r13:r13, r14:r14, r:rtsM[K].data[2].r, n:rtsM[K].data[2].n});

  // The middle position is left alone
  rtsW[K].push(rtsM[K].data[2]);

  theta = theta.neg();
  R = MP(I, MP(sM(theta.sin(), Q), sM(new Decimal(1).minus(theta.cos()), MM(Q,Q))));

  r12 = Mmul(R, rtsM[K].data[3].r12);
  r13 = Mmul(R, rtsM[K].data[3].r13);
  r14 = Mmul(R, rtsM[K].data[3].r14);

  rtsW[K].push({r12:r12, r13:r13, r14:r14, r:rtsM[K].data[2].r, n:rtsM[K].data[2].n});

  r12 = Mmul(R, Mmul(R, rtsM[K].data[4].r12));
  r13 = Mmul(R, Mmul(R, rtsM[K].data[4].r13));
  r14 = Mmul(R, Mmul(R, rtsM[K].data[4].r14));

  rtsW[K].push({r12:r12, r13:r13, r14:r14, r:rtsM[K].data[2].r, n:rtsM[K].data[2].n});

  return doTR(rtsW[K], DT);
}

function doTR(rts, DT)
{
  if (rts.data.length == 5)
  {
    let a12 = {
      x: rts.data[0].r12.x.neg().plus(rts.data[1].r12.x.times(16)).minus(rts.data[2].r12.x.times(30)).plus(rts.data[3].r12.x.times(16)).minus(rts.data[4].r12.x).div(DT.pow(2).times(12)),
      y: rts.data[0].r12.y.neg().plus(rts.data[1].r12.y.times(16)).minus(rts.data[2].r12.y.times(30)).plus(rts.data[3].r12.y.times(16)).minus(rts.data[4].r12.y).div(DT.pow(2).times(12)),
      z: rts.data[0].r12.z.neg().plus(rts.data[1].r12.z.times(16)).minus(rts.data[2].r12.z.times(30)).plus(rts.data[3].r12.z.times(16)).minus(rts.data[4].r12.z).div(DT.pow(2).times(12))
    };

    let a13 = {
      x: rts.data[0].r13.x.neg().plus(rts.data[1].r13.x.times(16)).minus(rts.data[2].r13.x.times(30)).plus(rts.data[3].r13.x.times(16)).minus(rts.data[4].r13.x).div(DT.pow(2).times(12)),
      y: rts.data[0].r13.y.neg().plus(rts.data[1].r13.y.times(16)).minus(rts.data[2].r13.y.times(30)).plus(rts.data[3].r13.y.times(16)).minus(rts.data[4].r13.y).div(DT.pow(2).times(12)),
      z: rts.data[0].r13.z.neg().plus(rts.data[1].r13.z.times(16)).minus(rts.data[2].r13.z.times(30)).plus(rts.data[3].r13.z.times(16)).minus(rts.data[4].r13.z).div(DT.pow(2).times(12))
    };

    let a14 = {
      x: rts.data[0].r14.x.neg().plus(rts.data[1].r14.x.times(16)).minus(rts.data[2].r14.x.times(30)).plus(rts.data[3].r14.x.times(16)).minus(rts.data[4].r14.x).div(DT.pow(2).times(12)),
      y: rts.data[0].r14.y.neg().plus(rts.data[1].r14.y.times(16)).minus(rts.data[2].r14.y.times(30)).plus(rts.data[3].r14.y.times(16)).minus(rts.data[4].r14.y).div(DT.pow(2).times(12)),
      z: rts.data[0].r14.z.neg().plus(rts.data[1].r14.z.times(16)).minus(rts.data[2].r14.z.times(30)).plus(rts.data[3].r14.z.times(16)).minus(rts.data[4].r14.z).div(DT.pow(2).times(12))
    };

    let r12=rts.data[2].r12;
    let r13=rts.data[2].r13;
    let r14=rts.data[2].r14;

    if (corr == 1 && ("n" in rts.data[2]))
    {
      let n = rts.data[2].n;
      let r = rts.data[2].r;

      let da12 = smul(GMD.times(-3).div(r.pow(4)), vadd(smul(norm(r12).pow(2).times(1.5).minus(dot(n, r12).pow(2).times(2.5)), n), X(r12, X(r12, n))));
      let da13 = smul(GMD.times(-3).div(r.pow(4)), vadd(smul(norm(r13).pow(2).times(1.5).minus(dot(n, r13).pow(2).times(2.5)), n), X(r13, X(r13, n))));
      let da14 = smul(GMD.times(-3).div(r.pow(4)), vadd(smul(norm(r14).pow(2).times(1.5).minus(dot(n, r14).pow(2).times(2.5)), n), X(r14, X(r14, n))));

      a12 = vadd(a12, da12);
      a13 = vadd(a13, da13);
      a14 = vadd(a14, da14);
    }

    let V2 = dot(r12,X(r13,r14));
    let V3 = dot(r13,X(r14,r12));
    let V4 = dot(r14,X(r12,r13));
    let V = V2.plus(V3).plus(V4).div(3);

    let tr = dot(a12,X(r13,r14)).plus(dot(a13,X(r14,r12))).plus(dot(a14,X(r12,r13))).div(V);

    return tr;
  }
  return new Decimal(NaN);
}

function eccentricity(state, refstate)
{
  let r = Decimal.sqrt(state.x.plus(refstate.x).pow(2).plus(state.y.plus(refstate.y).pow(2)).plus(state.z.plus(refstate.z).pow(2))).toNumber();
  let v = Decimal.sqrt(state.vx.pow(2).plus(state.vy.pow(2)).plus(state.vz.pow(2))).toNumber();
  let h = Math.sqrt(GM*r);
  let e = (v*v/GM - 1/r)/(v*v/GM + 1/r);
  return e;
}

function avgeccentricity(states, refstate)
{
  let x = new Decimal(0); let y = new Decimal(0); let z = new Decimal(0);
  let vx = new Decimal(0); let vy = new Decimal(0); let vz = new Decimal(0);
  for (let i=0; i<states.length; i++)
  {
    x = x.plus(states[i].x);
    y = y.plus(states[i].y);
    z = z.plus(states[i].z);
    vx = vx.plus(states[i].vx);
    vy = vy.plus(states[i].vy);
    vz = vz.plus(states[i].vz);
  }
  x = x.div(states.length).plus(refstate.x);
  y = y.div(states.length).plus(refstate.y);
  z = z.div(states.length).plus(refstate.z);
  vx = vx.div(states.length);
  vy = vy.div(states.length);
  vz = vz.div(states.length);

  let hx = y.times(vz).minus(z.times(vy)).toNumber();
  let hy = z.times(vx).minus(x.times(vz)).toNumber();
  let hz = x.times(vy).minus(y.times(vx)).toNumber();
  let h = Math.sqrt(hx*hx + hy*hy + hz*hz);
  let r = Decimal.sqrt(x.pow(2).plus(y.pow(2)).plus(z.pow(2))).toNumber();
  let v = Decimal.sqrt(vx.pow(2).plus(vy.pow(2)).plus(vz.pow(2))).toNumber();
  let e = Math.sqrt(1 + h*h/GM/GM*(v*v-2*GM/r));
  return e;
}

function orbitalPeriod(stateVector, refstate)
{
  let r = Decimal.sqrt(stateVector.x.plus(refstate.x).pow(2).plus(stateVector.y.plus(refstate.y).pow(2)).plus(stateVector.z.plus(refstate.z).pow(2))).toNumber();
  let v = Decimal.sqrt(stateVector.vx.pow(2).plus(stateVector.vy.pow(2)).plus(stateVector.vz.pow(2))).toNumber();
  let epsilon = 0.5*v*v - GM/r;
  let a = -GM/2/epsilon;
  let T = 2*Math.PI*Math.sqrt(a*a*a/GM);
  return T;
}

function meanOrbitalPeriod(stateVectors, refstate)
{
  let sumX = new Decimal(0);
  let sumY = new Decimal(0);
  let sumZ = new Decimal(0);
  let sumVx = new Decimal(0);
  let sumVy = new Decimal(0);
  let sumVz = new Decimal(0);
  for (let i=0; i<stateVectors.length; i++) {
    sumX = sumX.plus(stateVectors[i].x);
    sumY = sumY.plus(stateVectors[i].y);
    sumZ = sumZ.plus(stateVectors[i].z);
    sumVx = sumVx.plus(stateVectors[i].vx);
    sumVy = sumVy.plus(stateVectors[i].vy);
    sumVz = sumVz.plus(stateVectors[i].vz);
  }
  let mean = {
    x: sumX.div(stateVectors.length),
    y: sumY.div(stateVectors.length),
    z: sumZ.div(stateVectors.length),
    vx: sumVx.div(stateVectors.length),
    vy: sumVy.div(stateVectors.length),
    vz: sumVz.div(stateVectors.length)
  };
  return orbitalPeriod(mean, refstate);
}

var inPROP = false;

function propagate(dontMove = false, forward = true)
{
  if (inPROP) return;
  inPROP = true;

  const DT = forward ? new Decimal(step) : new Decimal(step).neg();
  let NT = Math.floor(86400/step);
  if (NT > 96) NT = 96;
  if (NT < 1) NT = 1;

  if (!dontMove)
  {
    for (let i = 0; i < NT; i++)
    {
      traceT = new Decimal(0);
      traceW = new Decimal(0);

      stdevT = new Decimal(0);
      stdevW = new Decimal(0);

      for (let j = 0; j < states.length; j++)
      {
        states[j] = integrator(states[j], DT, refstate);
      }
      var newstate = integrator(refstate, DT, nullstate);
      for (let j = 0; j < states.length; j++)
      {
        states[j].x = states[j].x.plus(refstate.x).minus(newstate.x);
        states[j].y = states[j].y.plus(refstate.y).minus(newstate.y);
        states[j].z = states[j].z.plus(refstate.z).minus(newstate.z);
      }
      refstate = newstate;

      for (let K = 0; K < 4; K++)
      {
        let tT = trT(states, DT, K);
        let tM = trM(states, DT, K); // Needed by trW
        let tW = trW(DT, K); // Uses rtsT and rtsM, created by trT and trM

        traceT = traceT.plus(tT);
        traceW = traceW.plus(tW);

        stdevT = stdevT.plus(tT.pow(2));
        stdevW = stdevW.plus(tW.pow(2));

        // We need to exclude outliers...

        const Z = new Decimal(5);

        if (isFinite(tW.toNumber()) && (Wnum < 10 || tW.abs().minus(Wsum.div(Wnum).abs()).lt(Z.times(Decimal.sqrt(Wdev.div(Wnum).minus(Wsum.div(Wnum).pow(2))))))) {
          Wsum = Wsum.plus(tW);
          Wdev = Wdev.plus(tW.pow(2));
          Wnum++;
        }
        if (isFinite(tT.toNumber()) && (Tnum < 10 || tT.abs().minus(Tsum.div(Tnum).abs()).lt(Z.times(Decimal.sqrt(Tdev.div(Tnum).minus(Tsum.div(Tnum).pow(2))))))) {
          Tsum = Tsum.plus(tT);
          Tdev = Tdev.plus(tT.pow(2));
          Tnum++;
        }
      }
      traceT = traceT.div(4);
      traceW = traceW.div(4);

      stdevT = Decimal.sqrt(stdevT.div(4).minus(traceT.pow(2)));
      stdevW = Decimal.sqrt(stdevW.div(4).minus(traceW.pow(2)));

      time += DT.toNumber() / 86400;
    }
  }
  let theVolume = volume(transformCoordinates(states, refstate)).toFixed(2);

  if (doCSV && rtsT[0].data.length > 2)
  {
    let ACOS = function(r1, r2)
    {
      return new Decimal(180).div(Decimal.acos(dot(r1, r2).div(norm(r1).times(norm(r2))))).toNumber();
    }

    let rMin = new Decimal("1e99");
    let rMax = new Decimal(0);
    let aMin = new Decimal(180);
    let aMax = new Decimal(0);
    for (let K = 0; K < 4; K++)
    {
      for (let prop of ['r12', 'r13', 'r14'])
      {
        let r = norm(rtsT[K].data[2][prop]);
        if (r.lt(rMin)) rMin = r;
        if (r.gt(rMax)) rMax = r;
      }
      for (let prop of [['r12','r13'],['r13','r14'],['r14','r12']])
      {
        let a = ACOS(rtsT[K].data[2][prop[0]], rtsT[K].data[2][prop[1]]);
        if (new Decimal(a).lt(aMin)) aMin = new Decimal(a);
        if (new Decimal(a).gt(aMax)) aMax = new Decimal(a);
      }
    }

    strLog += "" + time + ", " + theVolume + ", " + traceT.toNumber() + ", " +
               stdevT.toNumber() + ", " + traceW.toNumber() + ", " + stdevW.toNumber() + ", " +
               rMin.toNumber() + ", " + rMax.toNumber() + ", " + aMin.toNumber() + ", " + aMax.toNumber() + ", " +
               norm(rtsT[0].data[2].r12).toNumber() + ", " + norm(rtsT[0].data[2].r13).toNumber() + ", " + norm(rtsT[0].data[2].r14).toNumber() + ", " +
               norm(rtsT[1].data[2].r12).toNumber() + ", " + norm(rtsT[1].data[2].r13).toNumber() + ", " + norm(rtsT[2].data[2].r12).toNumber() +
               rtsT.reduce((a,b,i) => a + ", " + ACOS(b.data[2].r12,b.data[2].r13) + ", " + ACOS(b.data[2].r13,b.data[2].r14) + ", " + ACOS(b.data[2].r14,b.data[2].r12),"") + ", " + (rtsT[0].data[2].r.div(AU)).toNumber() + "\n";
  }

  if (theVolume[0] != '-') theVolume = "&nbsp;" + theVolume;
  let strTraceT = traceT.toPrecision(5);
  let strTraceW = traceW.toPrecision(5);
  let strStdevT = "&plusmn;" + stdevT.toPrecision(5);
  let strStdevW = "&plusmn;" + stdevW.toPrecision(5);

  trTavg = Tsum.div(Tnum).toPrecision(5);
  trWavg = Wsum.div(Wnum).toPrecision(5);
  trTdev = "&plusmn;" + Decimal.sqrt(Tdev.div(Tnum).minus(Tsum.div(Tnum).pow(2))).toPrecision(5);
  trWdev = "&plusmn;" + Decimal.sqrt(Wdev.div(Wnum).minus(Wsum.div(Wnum).pow(2))).toPrecision(5);

  if (strTraceT[0] != '-') strTraceT = "&nbsp;" + strTraceT;
  if (strTraceW[0] != '-') strTraceW = "&nbsp;" + strTraceW;

  if (trTavg[0] != '-') trTavg = "&nbsp;" + trTavg;
  if (trWavg[0] != '-') trWavg = "&nbsp;" + trWavg;

  document.getElementById("volume").innerHTML = theVolume;
  document.getElementById("trT").innerHTML = strTraceT + "<br/>" + strStdevT;
  document.getElementById("trW").innerHTML = strTraceW + "<br/>" + strStdevW;

  document.getElementById("trTavg").innerHTML = trTavg + "<br/>" + trTdev;
  document.getElementById("trWavg").innerHTML = trWavg + "<br/>" + trWdev;

  document.getElementById("time").innerText = parseFloat(time).toFixed(2);

  render();

  inPROP = false;
}

function drawShinySphere(ctx, x, y, radius, r, g, b)
{
  // Create a radial gradient
  const gradient = ctx.createRadialGradient(x, y, 0, x, y, radius);

  // Add color stops to create the shiny effect
  gradient.addColorStop(0, `rgba(${r}, ${g}, ${b}, 1)`); // Bright center
  gradient.addColorStop(0.8, `rgba(${Math.round(r/2.5)}, ${Math.round(g/2.5)}, ${Math.round(b/2.5)}, 0.8)`); // Slightly darker
  gradient.addColorStop(1, 'rgba(0, 0, 0, 0)'); // Faded edges

  // Set the fill style to the radial gradient
  ctx.fillStyle = gradient;

  // Draw the circle
  ctx.beginPath();
  ctx.arc(x, y, radius, 0, 2 * Math.PI);
  ctx.closePath();
  ctx.fill();
}

function setView(v)
{
  var v = document.getElementById("view").value;
  var d = 3 + 0.27 * document.getElementById("zoom").value;

  if (v == 1) { camera.x = d; camera.y = 0; camera.z = 0; }
  if (v == 2) { camera.x = 0; camera.y = d; camera.z = 0; }
  if (v == 3) { camera.x = 0; camera.y = 0; camera.z = d; }
  camera.phi = camera.theta = 0;
  view = v;
  if (!working) render();
}

function render()
{
  var spacecrafts = transformCoordinates(states, refstate);
  var v = document.getElementById("view").value * 1;

  spacecrafts.forEach(sc =>
  {
    var tempX, tempY, tempZ;

    switch (v)
    {
    case 1:
      tempX = sc.x.times(Decimal.cos(camera.phi)).plus(sc.y.times(Decimal.sin(camera.phi))).toNumber();
      sc.y = sc.x.times(Decimal.sin(camera.phi)).neg().plus(sc.y.times(Decimal.cos(camera.phi))).toNumber();
      sc.z = sc.z.times(Decimal.cos(camera.theta)).plus(new Decimal(tempX).times(Decimal.sin(camera.theta))).toNumber();
      sc.x = tempX;
      break;
    case 2:
      tempY = sc.y.times(Decimal.cos(camera.theta)).minus(sc.x.times(Decimal.sin(camera.theta))).toNumber();
      sc.x = sc.y.times(Decimal.sin(camera.theta)).plus(sc.x.times(Decimal.cos(camera.theta))).toNumber();
      sc.z = sc.z.times(Decimal.cos(camera.phi)).minus(new Decimal(tempY).times(Decimal.sin(camera.phi))).toNumber();
      sc.y = tempY;
      break;
    case 3:
      tempZ = sc.z.times(Decimal.cos(camera.theta)).minus(sc.y.times(Decimal.sin(camera.theta))).toNumber();
      sc.y = sc.z.times(Decimal.sin(camera.theta)).plus(sc.y.times(Decimal.cos(camera.theta))).toNumber();
      sc.x = sc.x.times(Decimal.cos(camera.phi)).minus(new Decimal(tempZ).times(Decimal.sin(camera.phi))).toNumber();
      sc.z = tempZ;
      break;
    }
  });

  // Clear the canvas
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  if (view == 3)
  {
    // Draw sun
    ctx.fillStyle = 'yellow';
    ctx.beginPath();
    let cw = Math.min(canvas.width, canvas.height);
    ctx.arc(canvas.width/2 + camera.phi*4/Math.PI*cw/2, canvas.height/2 + camera.theta*4/Math.PI*cw/2, -5*scale/sunZ.toNumber(), 0, 2 * Math.PI);
    ctx.fill();
  }

  // Sort the indices array based on the 'z' property of the corresponding objects in the original array

  if (view == 1) indices.sort((a, b) => spacecrafts[a].x - spacecrafts[b].x);
  if (view == 2) indices.sort((a, b) => spacecrafts[a].y - spacecrafts[b].y);
  if (view == 3) indices.sort((a, b) => spacecrafts[a].z - spacecrafts[b].z);

  var views = [];

  // Iterate through the sorted indices array and access the original array elements
  indices.forEach(index =>
  {
    var spacecraft = spacecrafts[index];

    // Calculate distance from the camera
    const dx = spacecraft.x - camera.x;
    const dy = spacecraft.y - camera.y;
    const dz = spacecraft.z - camera.z;

    const distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

    // Calculate sphere size based on distance
    const size = 500 / Math.pow(distance, 1.5);

    // Calculate screen coordinates
    var screenX, screenY;
    var tempX, tempY, tempZ;

    const canvSize = Math.min(canvas.width, canvas.height);

    if (view == 1)
    {
      if (dx > 0) return;
      tempX = (spacecraft.y / (camera.x - spacecraft.x)) * canvSize;
      tempY = (spacecraft.z / (camera.x - spacecraft.x)) * canvSize;
    }
    if (view == 2)
    {
      if (dy > 0) return;
      tempX = (spacecraft.z / (camera.y - spacecraft.y)) * canvSize;
      tempY = (spacecraft.x / (camera.y - spacecraft.y)) * canvSize;
    }
    if (view == 3)
    {
      if (dz > 0) return;
      tempX = (spacecraft.x / (camera.z - spacecraft.z)) * canvSize;
      tempY = (spacecraft.y / (camera.z - spacecraft.z)) * canvSize;
    }

    screenX = canvas.width / 2 + tempX;
    screenY = canvas.height / 2 - tempY;

    // Draw sphere
    views.push({X:screenX, Y:screenY, S:size, R:spacecraft.r, G:spacecraft.g, B:spacecraft.b});
  });

  // Draw lines connecting pairs of spacecraft
  for (let i = 0; i < views.length; i++)
  {
    const spacecraft1 = views[i];
    drawShinySphere(ctx, spacecraft1.X, spacecraft1.Y, spacecraft1.S, spacecraft1.R, spacecraft1.G, spacecraft1.B);
    for (let j = i + 1; j < views.length; j++)
    {
      const spacecraft2 = views[j];

      // Create a gradient for the line
      const gradient = ctx.createLinearGradient(spacecraft1.X, spacecraft1.Y, spacecraft2.X, spacecraft2.Y);
      gradient.addColorStop(0, `rgba(${spacecraft1.R}, ${spacecraft1.G}, ${spacecraft1.B}, ${Math.pow(spacecraft1.S,3) / 5000})`);
      gradient.addColorStop(1, `rgba(${spacecraft2.R}, ${spacecraft2.G}, ${spacecraft2.B}, ${Math.pow(spacecraft2.S,3) / 5000})`);

      // Draw the line
      ctx.beginPath();
      ctx.moveTo(spacecraft1.X, spacecraft1.Y);
      ctx.lineTo(spacecraft2.X, spacecraft2.Y);
      ctx.strokeStyle = gradient;
      ctx.lineWidth = 1;
      ctx.stroke();
    }
  }

  topView.update(states);

  document.getElementById('phi').innerText = new Decimal(-camera.phi).times(180).div(Math.PI).toFixed(2);
  document.getElementById('theta').innerText = new Decimal(camera.theta).times(180).div(Math.PI).toFixed(2);
}

function zoom()
{
  var d = 3 + 0.27 * document.getElementById("zoom").value;
  if (view == 1) camera.x = d;
  if (view == 2) camera.y = d;
  if (view == 3) camera.z = d;
  if (!working) render();
}

var worker = null;
var working = false;

function doWork(e)
{
  self.processing = false;
  if (e.data.cmd === 'run')
  {
    if (!self.hasOwnProperty('timer') || timer == 0)
      self.timer = setInterval(function ()
      {
        if (!self.processing)
        {
          self.processing = true;
          self.postMessage(0);
        }
      }, 50);
  }
  if (e.data.cmd == 'done')
  {
    self.processing = false;
  }
  if (e.data.cmd === 'stop')
  {
    if (self.hasOwnProperty('timer')) clearInterval(timer);
    self.timer = 0;
  }
}

function start()
{
  document.getElementById('skipForward').disabled = true;
  document.getElementById('skipBack').disabled = true;
  working = true;
  if (worker == null)
  {
    worker = new Worker(window.URL.createObjectURL(new Blob(["onmessage=" + doWork.toString()], {type: "text/javascript"})));
    worker.onmessage = function(e)
    {
      if (working) propagate();
      worker.postMessage({cmd: 'done'});
    }
  }
  worker.postMessage({cmd: 'run'});
}

function stop()
{
  document.getElementById('skipForward').disabled = false;
  document.getElementById('skipBack').disabled = false;
  if (worker) worker.postMessage({cmd: 'stop'});
  working = false;
}

function play()
{
  document.getElementById("init").disabled = true;
  if (!working) start();
  else stop();
}

function skipForward()
{
  document.getElementById("init").disabled = true;
  if (!working) propagate(false, true);
}

function skipBack()
{
  document.getElementById("init").disabled = true;
  if (!working) propagate(false, false);
}

function onProp()
{
  if (!working) propagate(true);
}

var canvas;
var ctx;

function onInit(init)
{
  if (init >= 0 && init < spacecraftsets.length)
  {
    spacecrafts = spacecraftsets[init].states;
    topView.bufferSize = spacecraftsets[init].tail;
    refstate = new State(spacecraftsets[init].refstate.x, spacecraftsets[init].refstate.y,
                         spacecraftsets[init].refstate.z, spacecraftsets[init].refstate.vx,
                         spacecraftsets[init].refstate.vy, spacecraftsets[init].refstate.vz, 1, 1, 1);
    scale = spacecraftsets[init].scale;
    step = spacecraftsets[init].step;
  }

  states = spacecrafts.map(spacecraft => new State(spacecraft.x, spacecraft.y, spacecraft.z,
                                                   spacecraft.vx, spacecraft.vy, spacecraft.vz,
                                                   spacecraft.r, spacecraft.g, spacecraft.b));

  // Create an array of indices
  indices = spacecrafts.map((_, index) => index);

  document.getElementById("e").innerText = avgeccentricity(states, refstate).toFixed(4);
  document.getElementById("T").innerText = (meanOrbitalPeriod(states, refstate)/86400).toFixed(2);
  propagate(true);
}

window.addEventListener('DOMContentLoaded', () =>
{
  const select = document.getElementById("init");
  for (let i = 0; i < spacecraftsets.length; i++)
  {
    const option = document.createElement('option');
    option.value = i;
    option.text = spacecraftsets[i].name;
    select.appendChild(option);
  }
  select.value = init;

  doResize();
  onInit(init);
  document.getElementById("view").value = view;
  setView(view);
});


function doResize()
{
  const body = document.querySelector('body');
  var ratioX = window.innerWidth / 1200.0;
  var ratioY = window.innerHeight / 800.0;
  var ratio = Math.min(ratioX, ratioY);
  if (ratio < 0.4) ratio = 0.4;
  if (ratio > 3) ratio = 3;
  body.style.zoom = ratio;

  const visualization = document.getElementById('visualization');
  visualization.innerHTML = "";
  canvas = document.createElement('canvas');
  ctx = canvas.getContext('2d');
  canvas.width = visualization.clientWidth;
  canvas.height = visualization.clientHeight;
  visualization.appendChild(canvas);

  let z = ([camera.x, camera.y, camera.z])[view - 1];
  document.getElementById("zoom").value=(z - 3)*100/27;

  function doZoom(delta, step=5)
  {
    var s = document.getElementById("zoom");
    var z = s.value;
    z = z - z%step;
    if (delta < 0 && z > 0) z -= step;
    if (delta > 0 && z < 100) z += step;
    s.value = z;
    zoom();
  }

  // Add a scroll event listener to the canvas
  canvas.addEventListener("wheel", function(event)
  {
    event.preventDefault();
    doZoom(event.deltaY, 5);
  });

  function doClick(e)
  {
    if (!canvas.wasdragging && !canvas.wasscaling) play();
    canvas.wasdragging = false;
    canvas.wasscaling = false;
  }

  function distance(touches)
  {
    return Math.sqrt((touches[1].clientX - touches[0].clientX)**2
                     + (touches[1].clientY - touches[0].clientY)**2);
  }

  canvas.addEventListener("click", doClick);

  function doStart(e)
  {
    canvas.dragging = true;
    if (e.touches)
    {
      canvas.lastX = e.touches[0].clientX;
      canvas.lastY = e.touches[0].clientY;
      if (e.touches.length === 2)
      {
        canvas.lastD = distance(e.touches);
      }
    }
    else
    {
    canvas.lastX = e.offsetX;
    canvas.lastY = e.offsetY;
    }
  }

  canvas.addEventListener("mousedown", doStart);
  canvas.addEventListener("touchstart", (e) =>
  {
    e.preventDefault();
    doStart(e);
  });

  function doMove(e)
  {
    if (canvas.dragging)
    {
      let dx = 0;
      let dy = 0;

      if (e.touches)
      {
        if (e.touches.length === 2)
        {
          let d = distance(e.touches);
          if (canvas.lastD != d)
          {
            doZoom(canvas.lastD - d, 1);
            canvas.wasscaling = true;
            canvas.lastD = d;
          }
        }
        else
        {
          dx = e.touches[0].clientX - canvas.lastX;
          dy = e.touches[0].clientY - canvas.lastY;

          canvas.lastX = e.touches[0].clientX;
          canvas.lastY = e.touches[0].clientY;
          if (dx != 0 || dy != 0) canvas.wasdragging = true;
        }
      }
      else
      {
        dx = e.offsetX - canvas.lastX;
        dy = e.offsetY - canvas.lastY;

        canvas.lastX = e.offsetX;
        canvas.lastY = e.offsetY;
        canvas.wasdragging = true;
      }

      camera.phi -= .5 * dx / 180 * Math.PI;
      camera.theta -= .5 * dy / 180 * Math.PI;

      if (camera.phi < -0.5*Math.PI) camera.phi = -0.5*Math.PI;
      if (camera.phi > 0.5*Math.PI) camera.phi = 0.5*Math.PI;
      if (camera.theta < -0.5*Math.PI) camera.theta = -0.5*Math.PI;
      if (camera.theta > 0.5*Math.PI) camera.theta = 0.5*Math.PI;
      if (!working) render();
    }
  }

  canvas.addEventListener("mousemove", doMove);
  canvas.addEventListener("touchmove", (e) =>
  {
    e.preventDefault();
    doMove(e);
  });

  function doStop(e)
  {
    canvas.dragging = false;
  }

  canvas.addEventListener("mouseup", doStop);
  canvas.addEventListener("touchend", (e) =>
  {
    doStop();
    doClick();
  });

  topView.init(canvas, canvas.width-topView.width, canvas.height-topView.height);
}

let topView =
{
  width: 200,
  height: 200,
  context: null,
  buffer: [],
  bufferSize: 300,
  x: 0, // x coordinate of top left of view
  y: 0, // y coordinate of top left of view

  init(canvas, x, y)
  {
    this.context = canvas.getContext('2d');
    this.x = x;
    this.y = y;
  },

  update(stateVectors)
  {
    // Average state vectors
    let avgX = 0;
    let avgY = 0;
    let avgZ = 0;
    stateVectors.forEach(v =>
    {
      avgX += v.x.toNumber();
      avgY += v.y.toNumber();
      avgZ += v.z.toNumber();
    });
    avgX = avgX / stateVectors.length + refstate.x.toNumber();
    avgY = avgY / stateVectors.length + refstate.y.toNumber();
    avgZ = avgZ / stateVectors.length + refstate.z.toNumber();
    let avgState = {x: avgX / AU, y: -avgY / AU, z: avgZ / AU};

    // Add latest state to buffer
    if (working)
    {
      this.buffer.push(avgState);
      if (this.buffer.length > this.bufferSize)
      {
        this.buffer.shift();
      }
    }

    // Clear view area only
    this.context.clearRect(this.x, this.y, this.width, this.height);

    // Draw sun
    this.context.fillStyle = 'yellow';
    this.context.beginPath();
    this.context.arc(this.x+this.width/2, this.y+this.height/2, 8, 0, 2 * Math.PI);
    this.context.fill();

    // Scale state vector coordinates to view
    let x = (avgState.x * AU / scale + 2.5) * this.width / 5 + this.x;
    let y = (avgState.y * AU / scale + 2.5) * this.height / 5 + this.y;

    // Draw current position
    this.context.fillStyle = '#8FF';
    this.context.beginPath();
    this.context.arc(x, y, 4, 0, 2 * Math.PI);
    this.context.fill();

    // Draw trail
    this.context.lineWidth = 1;
    this.context.beginPath();
    this.context.moveTo(x, y);
    let tstep = Math.floor(this.buffer.length / 100);
    if (tstep < 1) tstep = 1;
    for (let i = this.buffer.length - 2; i >= 0; i -= tstep)
    {
      let state = this.buffer[i];
      let x = (state.x * AU / scale + 2.5) * this.width / 5 + this.x;
      let y = (state.y * AU / scale + 2.5) * this.height / 5 + this.y;
      this.context.lineTo(x, y);
      // Fade color
      this.context.strokeStyle = `rgba(224, 224, 255, ${(i/this.buffer.length)**0.125})`;
      this.context.stroke();
      this.context.beginPath();
      this.context.moveTo(x, y);
    }

    // Draw border around view area
    this.context.strokeStyle = 'lightgray';
    this.context.strokeRect(this.x, this.y, this.width, this.height);

    this.context.stroke();
  }
}

window.addEventListener('resize', () => { doResize(); onProp(); });
