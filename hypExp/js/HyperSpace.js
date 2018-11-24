//-------------------------------------------------------
// Global Variables
//-------------------------------------------------------
var g_renderer;
var g_effect;
var g_controls;
var g_geometry;
var g_currentBoost;
var g_stereoBoosts = [];
var g_leftCurrentBoost;
var g_rightCurrentBoost;
var g_cellBoost;
var g_invCellBoost;
var g_screenResolution;
var g_controllerBoosts = [];
var g_controllerDualPoints = [];
var g_vr;
var g_composer;

//-------------------------------------------------------
// Scene Variables
//-------------------------------------------------------
var maxSteps = 50;
var textFPS;
var time;

//-------------------------------------------------------
// FPS Manager
//-------------------------------------------------------
var m_stepDamping = 0.75;
var m_stepAccum = 0;
var fpsLog = new Array(10);
fpsLog.fill(g_targetFPS.value);

var fps = {
	lastTime: null,
	getFPS: function () {
		if(!this.lastTime) {
			this.lastTime = new Date();
			return null;
		}
		var date = new Date();
		var currentFps = 1000 / (date - this.lastTime);
		this.lastTime = date;
		return currentFps;
	}
}

var calcMaxSteps = function(lastFPS, lastMaxSteps)
{
  if(guiInfo.autoSteps){
	  if(!lastFPS)
		  return lastMaxSteps;

	  fpsLog.shift();
	  fpsLog.push(lastFPS);
	  var averageFPS = Math.average(fpsLog);
	  textFPS.innerHTML = averageFPS.toPrecision(3);

	  // We don't want the adjustment to happen too quickly (changing maxSteps every frame is quick!),
	  // so we'll let fractional amounts m_stepAccumulate until they reach an integer value.
	  var newVal = Math.pow((averageFPS / g_targetFPS.value), (1 / 20)) * lastMaxSteps;
	  var diff = newVal - lastMaxSteps;
	  if(Math.abs( m_stepAccum ) < 1)
	  {
		  m_stepAccum += diff;
		  m_stepAccum *= m_stepDamping;
		  return lastMaxSteps;
	  }

	  newVal = lastMaxSteps + m_stepAccum;
	  newVal = Math.round(Math.clamp(newVal, 31, 127));
	  m_stepAccum = 0;
	  return newVal;
  }
  else {
    return guiInfo.maxSteps;
  }
}

//-------------------------------------------------------
// Sets up precalculated values
//-------------------------------------------------------
var hCWH = 0.6584789485;
var hCWK = 0.5773502692;
var gens;
var invGens;
var hCDP = [];
var simplexMirrors = [];

var initGenerators = function( p, q, r ){
  g_geometry = GetGeometry( p, q, r );
  var isCubical = p == 4 && q == 3;

  if( isCubical )
  {
    var invHCWK = 1.0/hCWK;
    
    hCDP[0] = new THREE.Vector4(invHCWK,0.0,0.0,1.0);
    hCDP[1] = new THREE.Vector4(0.0,invHCWK,0.0,1.0);
    hCDP[2] = new THREE.Vector4(0.0,0.0,invHCWK,1.0);
    if( g_geometry != Geometry.Euclidean ) {
      for( var i=0; i<3; i++ )
        hCDP[i].geometryNormalize(g_geometry);
    }

    gens = createCubeGenerators(g_geometry);
    invGens = invCubeGenerators(gens);

    simplexMirrors = [];
    for(var i = 0; i<4; i++){
      simplexMirrors.push(new THREE.Vector4());
    }
  }
  else
  {
    simplexMirrors = SimplexFacetsKlein( p, q, r );
    invGens = SimplexInverseGenerators( p, q, r );

    // invGens needs to be length-6;
    for(var i = 0; i<2; i++){
      invGens.push(translateByVector(g_geometry, new THREE.Vector3(0.0,0.0,0.0)));
    }
  }

  for(var i = 0; i<6; i++){
    g_controllerDualPoints.push(new THREE.Vector4());
  }
}

var createCubeGenerators = function(g){
  var gen0 = translateByVector(g, new THREE.Vector3(2.0*hCWH,0.0,0.0));
  var gen1 = translateByVector(g, new THREE.Vector3(-2.0*hCWH,0.0,0.0));
  var gen2 = translateByVector(g, new THREE.Vector3(0.0,2.0*hCWH,0.0));
  var gen3 = translateByVector(g, new THREE.Vector3(0.0,-2.0*hCWH,0.0));
  var gen4 = translateByVector(g, new THREE.Vector3(0.0,0.0,2.0*hCWH));
  var gen5 = translateByVector(g, new THREE.Vector3(0.0,0.0,-2.0*hCWH));
  return [gen0, gen1, gen2, gen3, gen4, gen5];
}

var invCubeGenerators = function(genArr){
  return [genArr[1],genArr[0],genArr[3],genArr[2],genArr[5],genArr[4]];
}


//-------------------------------------------------------
// Sets up the lights
//-------------------------------------------------------
var lightPositions = [];
var lightIntensities = [];
var attnModel = 1;

var initLights = function(g){
  lightPositions = [];
  lightIntensities = [];
  PointLightObject(g, new THREE.Vector3(0,0,1), new THREE.Vector4(0,0,1,1));
  PointLightObject(g, new THREE.Vector3(1.2,0,0), new THREE.Vector4(1,0,0,1));
  PointLightObject(g, new THREE.Vector3(0,1.1,0), new THREE.Vector4(0,1,0,1));
  PointLightObject(g, new THREE.Vector3(-1,-1,-1), new THREE.Vector4(1,1,1,1));
  //Add light info for controllers
  lightIntensities.push(new THREE.Vector4(0.49, 0.28, 1.0, 2));
  lightIntensities.push(new THREE.Vector4(1.0, 0.404, 0.19, 2));
}

//-------------------------------------------------------
// Sets up global objects
//-------------------------------------------------------
var globalObjectBoosts = [];
var invGlobalObjectBoosts = [];
var globalObjectRadii = [];
var globalObjectTypes = [];

//TODO: CREATE GLOBAL OBJECT CONSTRUCTORS
var initObjects = function(g){
  globalObjectBoosts = [];
  invGlobalObjectBoosts = [];
  globalObjectRadii = [];
  globalObjectTypes = [];
  SphereObject(g, new THREE.Vector3(0.5,0,0), 0.2); // geometry, position, radius/radii
  EllipsoidObject(g, new THREE.Vector3(-0.5,0,0), new THREE.Vector3(1.0,0.7,0.5)); //radii must be less than one!
  for(var i = 2; i<4; i++){ // We need to fill out our arrays with empty objects for glsl to be happy
    EmptyObject();
  }
}

//-------------------------------------------------------
// Sets up  our raymarch shader pass
//-------------------------------------------------------
var raymarchPass = function(){
  var pass = new THREE.ShaderPass(THREE.hyper);
  
  //Our massive list of uniforms
  pass.uniforms.isStereo.value = g_vr;
  pass.uniforms.geometry.value = 3;
  pass.uniforms.screenResolution.value = g_screenResolution;
  pass.uniforms.fov.value = guiInfo.fov;
  pass.uniforms.invGenerators.value = invGens;
  pass.uniforms.currentBoost.value = g_currentBoost;
  pass.uniforms.stereoBoosts.value = g_stereoBoosts; // need to make array with leftCurrentBoost and right
  pass.uniforms.cellBoost.value = g_cellBoost;
  pass.uniforms.invCellBoost.value = g_invCellBoost;
  pass.uniforms.maxSteps.value = maxSteps;
  pass.uniforms.lightPositions.value = lightPositions;
  pass.uniforms.lightIntensities.value = lightIntensities;
  pass.uniforms.attnModel.value = attnModel;
  pass.uniforms.texture.value = new THREE.TextureLoader().load("images/concrete2.png");
  pass.uniforms.controllerCount.value = 0;
  pass.uniforms.controllerBoosts.value = g_controllerBoosts;
  pass.uniforms.globalObjectBoosts.value = globalObjectBoosts;
  pass.uniforms.invGlobalObjectBoosts.value = invGlobalObjectBoosts;
  pass.uniforms.globalObjectRadii.value = globalObjectRadii;
  pass.uniforms.globalObjectTypes.value = globalObjectTypes;
  pass.uniforms.halfCubeDualPoints.value = hCDP;
  pass.uniforms.halfCubeWidthKlein.value = hCWK;
  pass.uniforms.cut4.value = g_cut4;
  pass.uniforms.sphereRad.value = g_sphereRad;
  pass.uniforms.tubeRad.value = g_tubeRad;
  pass.uniforms.vertexPosition.value = g_vertexPosition;
  pass.uniforms.vertexSurfaceOffset.value = g_vertexSurfaceOffset;
  pass.uniforms.useSimplex.value = false;
  pass.uniforms.simplexMirrorsKlein.value = simplexMirrors;

  //Our list of defines
  //currently set as constant values is in the shader for easier debugging
 // pass.defines.NUM_LIGHTS.value = lightPositions.length;
 // pass.defines.NUM_OBJECTS = globalObjectBoosts.length;
  return pass;
}


//-------------------------------------------------------
// Sets up the scene
//-------------------------------------------------------
var init = function(){
  //Setup our THREE scene--------------------------------
	time = Date.now();
	textFPS = document.getElementById('fps');
  g_renderer = new THREE.WebGLRenderer();

  g_screenResolution = new THREE.Vector2(window.innerWidth, window.innerHeight);
  g_renderer.setSize(g_screenResolution.x, g_screenResolution.y);
  document.body.appendChild(g_renderer.domElement);

  g_controls = new THREE.Controls();
  g_effect = new THREE.VREffect();
  g_controllerBoosts.push(new THREE.Matrix4()); g_controllerBoosts.push(new THREE.Matrix4());
  g_currentBoost = new THREE.Matrix4(); // boost for camera relative to central cell
  g_cellBoost = new THREE.Matrix4(); // boost for the cell that we are in relative to where we started
  g_invCellBoost = new THREE.Matrix4();
  g_geometry = Geometry.Hyperbolic; // we start off hyperbolic
  initGenerators(4,3,6);
  initLights(g_geometry);
  initObjects(g_geometry);
  //Setup dat GUI --- UI.js
  initGui();

  //-------------------------------------------------------
  // "Post" Processing - Since we are not using meshes we actually 
  //                     don't need to do traditional rendering we 
  //                     can just use post processed effects
  //-------------------------------------------------------
  
  //Composer 
  //**********************************************
  g_composer = new THREE.EffectComposer(g_renderer);
  
  //Shader Passes 
  //**********************************************
  //Raymarch
  g_raymarch = raymarchPass();
  g_composer.addPass(g_raymarch);
  g_raymarch.renderToScreen = true;

  //Generator for controllerScaleMatrix on the glsl side
  //console.log(translateByVector(g_geometry, new THREE.Vector3(0,0,0.2)).multiply(scaleMatrix));
  
  animate();
}

//-------------------------------------------------------
// Where our scene actually renders out to screen
//-------------------------------------------------------
var animate = function(){
  requestAnimationFrame(animate);
  g_controls.update();
  maxSteps = calcMaxSteps(fps.getFPS(), maxSteps);
  //THREE.VRController.update();
  g_raymarch.uniforms.maxSteps.value = maxSteps;
  g_raymarch.uniforms.controllerCount.value = THREE.VRController.controllers.length;
  g_composer.render();
}

//-------------------------------------------------------
// Where the magic happens
//-------------------------------------------------------
if(mobileCheck()){
  window.location.replace("http://www.michaelwoodard.net/hypVR-Ray_m/")
}
else{
  init();
}