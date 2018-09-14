//-------------------------------------------------------
// Global Variables
//-------------------------------------------------------
var g_material;
var g_controls;
var g_currentBoost;
var g_cellBoost;
var g_invCellBoost;

//-------------------------------------------------------
// Scene Variables
//-------------------------------------------------------
var scene, camera, renderer, composer;
var stats;
//-------------------------------------------------------
// Sets up precalculated values
//-------------------------------------------------------
var hCWH = 0.6584789485;
var gens;
var invGens;

var createGenerators = function(){
  var gen0 = translateByVector(new THREE.Vector3(2.0*hCWH,0.0,0.0));
  var gen1 = translateByVector(new THREE.Vector3(-2.0*hCWH,0.0,0.0));
  var gen2 = translateByVector(new THREE.Vector3(0.0,2.0*hCWH,0.0));
  var gen3 = translateByVector(new THREE.Vector3(0.0,-2.0*hCWH,0.0));
  var gen4 = translateByVector(new THREE.Vector3(0.0,0.0,2.0*hCWH));
  var gen5 = translateByVector(new THREE.Vector3(0.0,0.0,-2.0*hCWH));
  return [gen0, gen1, gen2, gen3, gen4, gen5];
}

var invGenerators = function(genArr){
  return [genArr[1],genArr[0],genArr[3],genArr[2],genArr[5],genArr[4]];
}


//-------------------------------------------------------
// Sets up the global objects
//-------------------------------------------------------
var lightPositions = [];
var lightIntensities = [];
var globalObjectBoost;

var initObjects = function(){
  PointLightObject(new THREE.Vector3(1.2,0,0), new THREE.Vector4(1,0,0,1));
  PointLightObject(new THREE.Vector3(0,1.2,0), new THREE.Vector4(0,1,0,1));
  PointLightObject(new THREE.Vector3(0,0,1.2), new THREE.Vector4(0,0,1,1));
  PointLightObject(new THREE.Vector3(-1,-1,-1), new THREE.Vector4(1,1,1,1));
  globalObjectBoost = new THREE.Matrix4().multiply(translateByVector(new THREE.Vector3(-0.5,0,0)));
}

//-------------------------------------------------------
// Set up shader
//-------------------------------------------------------

var initShader = function(screenRes){
  var pass = new THREE.ShaderPass(THREE.RaymarchShader);
  pass.uniforms.screenResolution.value = screenRes;
  pass.uniforms.invGenerators.value = invGens;
  pass.uniforms.currentBoost.value = g_currentBoost;
  pass.uniforms.cellBoost.value = g_cellBoost;
  pass.uniforms.invCellBoost.value = g_invCellBoost;
  pass.uniforms.lightPositions.value = lightPositions;
  pass.uniforms.lightIntensities.value = lightIntensities;
  pass.uniforms.globalObjectBoost.value = globalObjectBoost;
  return pass;
}

//-------------------------------------------------------
// Sets up the scene
//-------------------------------------------------------
var init = function(){
  //Setup our THREE scene--------------------------------
  scene = new THREE.Scene();
  renderer = new THREE.WebGLRenderer();
  camera = new THREE.OrthographicCamera(-1,1,1,-1,1/Math.pow(2,53),1);

  var screenRes = new THREE.Vector2(window.innerWidth, window.innerHeight);
  renderer.setSize(screenRes.x, screenRes.y);
  document.body.appendChild(renderer.domElement);

  //Initialize varirables, objects, and stats
  stats = new Stats(); stats.showPanel(1); stats.showPanel(2); stats.showPanel(0); document.body.appendChild(stats.dom);
  g_controls = new THREE.Controls(); g_currentBoost = new THREE.Matrix4();  g_cellBoost = new THREE.Matrix4(); g_invCellBoost = new THREE.Matrix4();
  gens = createGenerators(); invGens = invGenerators(gens); initObjects();
  //Composer
  composer = new THREE.EffectComposer(renderer);
  //Render Passes
  //var renderPass = new THREE.RenderPass(scene, camera);
  //composer.addPass(renderPass);
  //Shader Passes
  var pass1 = initShader(screenRes);
  composer.addPass(pass1);
  pass1.renderToScreen = true;
  //Let's get rendering
  animate();
}

//-------------------------------------------------------
// Where our scene actually renders out to screen
//-------------------------------------------------------
var animate = function(){
  stats.begin();
  requestAnimationFrame(animate);
  composer.render();
  g_controls.update();
  stats.end();
}

//-------------------------------------------------------
// Where the magic happens
//-------------------------------------------------------
init();