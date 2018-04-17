var scene;
var renderer;
var effect;
var camera;
var virtCamera;
var mesh;
var geom;
var material;
var controls;
var currentBoost;

//-------------------------------------------------------
// Scene Manipulator Functions & Variables
//-------------------------------------------------------
var gens;
var invGens;
var maxSteps = 31;
var hCWH = 0.6584789485;
var hCWK = 0.5773502692;
var sphereRad = 1.0;
var horosphereSize = 2.6;

var createGenerators = function(){
  var gen0 = translateByVector(new THREE.Vector3( 2.0*hCWH, 0.0, 0.0));
  var gen1 = translateByVector(new THREE.Vector3(-2.0*hCWH, 0.0, 0.0));
  var gen2 = translateByVector(new THREE.Vector3(0.0,  2.0*hCWH, 0.0));
  var gen3 = translateByVector(new THREE.Vector3(0.0, -2.0*hCWH, 0.0));
  var gen4 = translateByVector(new THREE.Vector3(0.0, 0.0,  2.0*hCWH));
  var gen5 = translateByVector(new THREE.Vector3(0.0, 0.0, -2.0*hCWH));
  return [gen0, gen1, gen2, gen3, gen4, gen5];
}

var invGenerators = function(genArr){
  return [genArr[1],genArr[0],genArr[3],genArr[2],genArr[5],genArr[4]];
}

var fps = {
  lastTime: new Date().getTime(),
  //frameNum: 0,
  getFPS: function(){
    //this.frameNum++;
    var date = new Date().getTime();
    var deltaTime = (date-this.lastTime)/1000;
    this.lastTime = date;
    //var res = this.frameNum/deltaTime;
    /*if(deltaTime>1){
      this.start = new Date().getTime();
      this.frameNum=0;
    }*/
    return 1/deltaTime;
  }
}
var fpsLog = new Array(10);
fpsLog.fill(30.0);

var calcMaxSteps = function(targetFPS, lastFPS, lastMaxSteps){
  fpsLog.shift();
  fpsLog.push(lastFPS);

  var averageFPS = 0.0;
  for(var i=0; i<fpsLog.length; i++){
    averageFPS += fpsLog[i];
  }
  averageFPS /= fpsLog.length;
  //console.log(Math.floor(averageFPS));
  return Math.max(Math.min(Math.round(Math.pow((averageFPS/targetFPS),(1/10)) * lastMaxSteps),127),31);
}

//-------------------------------------------------------
// Sets up the scene
//-------------------------------------------------------
var init = function(){
  //Setup our THREE scene--------------------------------
  scene = new THREE.Scene();
  renderer = new THREE.WebGLRenderer();
  document.body.appendChild(renderer.domElement);
  effect = new THREE.VREffect(renderer);
  effect.setSize(window.innerWidth, window.innerHeight);
  camera = new THREE.OrthographicCamera(-1,1,1,-1,1/Math.pow(2,53),1);
  virtCamera = new THREE.PerspectiveCamera(60,1,0.1,1);
  virtCamera.position.z = 0.1;
  cameraOffset = new THREE.Vector3();
  controls = new THREE.VRControls(virtCamera);
  gens = createGenerators(hCWH);
  invGens = invGenerators(gens);
  currentBoost = new THREE.Matrix4(); // boost for camera relative to central cell
  cellBoost = new THREE.Matrix4(); // boost for the cell that we are in relative to where we started
  lightSourcePosition = new THREE.Vector4(0.0,0.0,0.9801960588,1.400280084); // position on hyperboloid of light source, is lorentzNormalize(0,0,.7,1)
  //Setup our material----------------------------------
  loadShaders();
}

var loadShaders = function(){ //Since our shader is made up of strings we can construct it from parts
  var loader = new THREE.FileLoader();
  loader.setResponseType('text')
  loader.load('shaders/fragment.glsl',function(main){
    loader.load('shaders/hyperbolicScene.glsl', function(scene){
      loader.load('shaders/hyperbolicMath.glsl', function(math){
        loader.load('shaders/globalsInclude.glsl', function(globals){
          finishInit(globals.concat(math).concat(scene).concat(main));
        });
      });
    });
  });
}

var finishInit = function(fShader){
//  console.log(fShader);
  material = new THREE.ShaderMaterial({
    uniforms:{
      screenResolution:{type:"v2", value:new THREE.Vector2(window.innerWidth, window.innerHeight)},
      cameraPos:{type:"v3", value:virtCamera.position},
      cameraQuat:{type:"v4", value:virtCamera.quaternion},
      fov:{type:"f", value:virtCamera.fov},
      generators:{type:"m4v", value:gens},
      invGenerators:{type:"m4v", value:invGens},
      currentBoost:{type:"m4", value:currentBoost},
      cellBoost:{type:"m4", value:cellBoost},
      lightSourcePosition:{type:"v4", value:lightSourcePosition},
      maxSteps:{type:"i", value:maxSteps},
      sceneIndex:{type:"i", value: 1},
      halfCubeWidthKlein:{type:"f", value: hCWK},
      sphereRad:{type:"f", value:sphereRad},
      horosphereSize:{type:"f", value:horosphereSize}
    },
    vertexShader: document.getElementById('vertexShader').textContent,
    fragmentShader: fShader,
    transparent:true
  });
  //Setup dat GUI
  var gui = new dat.GUI();
  gui.add(material.uniforms.sceneIndex, 'value',{Sphere_horosphere: 1, Sphere_plane: 2, Medial_surface: 3, Cube_planes: 4}).name("Scene");
  //Setup a "quad" to render on-------------------------
  geom = new THREE.BufferGeometry();
  var vertices = new Float32Array([
    -1.0, -1.0, 0.0,
     1.0, -1.0, 0.0,
     1.0,  1.0, 0.0,

    -1.0, -1.0, 0.0,
     1.0,  1.0, 0.0,
    -1.0,  1.0, 0.0
  ]);
  geom.addAttribute('position',new THREE.BufferAttribute(vertices,3));
  mesh = new THREE.Mesh(geom, material);
  scene.add(mesh);
  animate();
}

//-------------------------------------------------------
// Where our scene actually renders out to screen
//-------------------------------------------------------
var animate = function(){
  controls.update();
  maxSteps = calcMaxSteps(30, fps.getFPS(), maxSteps);
  material.uniforms.maxSteps.value = maxSteps;
  effect.render(scene, camera, animate);
}

//-------------------------------------------------------
// Where the magic happens
//-------------------------------------------------------
init();


//-------------------------------------------------------
// Event listeners
//-------------------------------------------------------
var onResize = function(){
  effect.setSize(window.innerWidth, window.innerHeight);
  if(material != null){
    material.uniforms.screenResolution.value.x = window.innerWidth;
    material.uniforms.screenResolution.value.y = window.innerHeight;
  }
}
window.addEventListener('resize', onResize, false);