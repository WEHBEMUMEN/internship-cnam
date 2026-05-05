/**
 * DEIM Benchmark — Application Core
 * Handles state, Three.js initialization, and coordinate mapping.
 */

const JET = [{t:0,r:.23,g:.51,b:.96},{t:.2,r:.02,g:.71,b:.83},{t:.4,r:.06,g:.73,b:.51},{t:.6,r:.98,g:.8,b:.08},{t:.8,r:.98,g:.45,b:.09},{t:1,r:.94,g:.27,b:.27}];
function jet(t){t=Math.max(0,Math.min(1,t));for(let i=0;i<JET.length-1;i++){const a=JET[i],b=JET[i+1];if(t>=a.t&&t<=b.t){const f=(t-a.t)/(b.t-a.t);return[a.r+f*(b.r-a.r),a.g+f*(b.g-a.g),a.b+f*(b.b-a.b)];}}return[.94,.27,.27];}

class Sparkline{constructor(c){this.c=c;this.ctx=c.getContext('2d');this.data=[];}
update(d){this.data=d;this.draw();}
draw(){const{c,ctx,data}=this;const W=c.width,H=c.height;ctx.clearRect(0,0,W,H);ctx.fillStyle='rgba(241,245,249,0.8)';ctx.fillRect(0,0,W,H);if(!data.length)return;const vals=data.map(d=>Math.max(d.norm,1e-14));const mn=Math.log10(Math.min(...vals)),mx=Math.log10(Math.max(...vals))+.5;const toY=v=>H-((Math.log10(Math.max(v,1e-14))-mn)/(mx-mn))*H;const step=W/Math.max(data.length-1,1);ctx.strokeStyle='rgba(0,0,0,0.05)';ctx.lineWidth=1;for(let e=Math.floor(mn);e<=Math.ceil(mx);e++){const y=toY(Math.pow(10,e));if(y<0||y>H)continue;ctx.beginPath();ctx.moveTo(0,y);ctx.lineTo(W,y);ctx.stroke();ctx.fillStyle='#64748b';ctx.font='8px monospace';ctx.fillText(`10^${e}`,2,y-2);}ctx.strokeStyle='#8b5cf6';ctx.lineWidth=1.5;ctx.lineJoin='round';ctx.beginPath();data.forEach((d,i)=>{const x=i*step,y=toY(d.norm);i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);});ctx.stroke();data.forEach((d,i)=>{ctx.fillStyle='#8b5cf6';ctx.beginPath();ctx.arc(i*step,toY(d.norm),2,0,Math.PI*2);ctx.fill();});if(data.length){const last=data[data.length-1].norm;ctx.fillStyle=last<1e-4?'#10b981':'#f59e0b';ctx.font='bold 9px monospace';ctx.fillText(`Res: ${last.toExponential(2)}`,4,H-6);}}}

const METHODS = {
    fom:      {label:'FOM',      color:'#334155'},
    galerkin: {label:'Galerkin', color:'#3b82f6'},
    ecsw:     {label:'ECSW',     color:'#10b981'}
};

class DEIMBenchmarkApp {
    constructor() {
        this.engine = new NURBS2D();
        this.solverFOM = new IGANonlinearSolver(this.engine);
        this.romEngine = new ROMEngine(this.solverFOM);
        this.ecswEngine = new ECSWEngine();

        this.patch = null;
        this.method = 'fom';
        this.loadMag = 200;
        this.k = 25;
        this.deimM = 40;
        this.meshLevel = 1;
        this.loadType = 'tip';
        this.showDOFs = true;
        this.isTrained = false;
        this.lastResult = null;
        this.lastFomTime = null;
        this.lastFomResult = null;
        this._timer = null;

        // DEIM Explorer State
        this.explorerStep = 1;
        this.isExplorerActive = false;
        this.sensorSpheres = new THREE.Group();
        this.activeElementMesh = null;

        this.scene = new THREE.Scene();
        this.loadArrows = new THREE.Group();
        this.scene.add(this.loadArrows);
        this.dofPoints = null;
        this.camera = new THREE.PerspectiveCamera(45, window.innerWidth/window.innerHeight, 0.1, 1000);
        this.renderer = new THREE.WebGLRenderer({antialias:true});
        this.deformedMesh = null;
        this.ghostMesh = null;

        this.initThree();
        this.initUI();
        this.initCharts();
        this.loadBenchmark();
    }

    initThree() {
        this.scene.background = new THREE.Color(0xffffff);
        this.camera.position.set(5, 1, 14);
        this.camera.lookAt(5, 1, 0);
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
        document.getElementById('canvas-container').appendChild(this.renderer.domElement);
        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableRotate = false;
        this.controls.screenSpacePanning = true;
        this.controls.mouseButtons = {
            LEFT: THREE.MOUSE.PAN,
            MIDDLE: THREE.MOUSE.DOLLY,
            RIGHT: THREE.MOUSE.ROTATE
        };
        this.controls.target.set(5, 1, 0);
        this.controls.update();
        this.controls.addEventListener('change', () => this._render());
        this.scene.add(new THREE.AmbientLight(0xffffff, 1.2));
        const dl = new THREE.DirectionalLight(0xffffff, 0.6);
        dl.position.set(5, 5, 10);
        this.scene.add(dl);
        window.addEventListener('resize', () => {
            this.camera.aspect = window.innerWidth/window.innerHeight;
            this.camera.updateProjectionMatrix();
            this.renderer.setSize(window.innerWidth, window.innerHeight);
            this._render();
        });
    }

    _render() { this.renderer.render(this.scene, this.camera); }

    loadBenchmark() {
        this.clearMesh();
        this.patch = NURBSPresets.generateCantilever(10, 2);
        RefineUtils.apply(this.engine, this.patch, { p: 3, h: this.meshLevel });
        this.solverFOM.E = 100000;
        this.solverFOM.nu = 0.3;
        this.solverFOM.thickness = 1.0;

        this.isTrained = false;
        this.method = 'fom';
        this.lastFomResult = null;
        this.lastFomTime = null;
        document.querySelectorAll('[data-method]').forEach(b => {
            b.classList.remove('active');
            if(b.dataset.method === 'fom') b.classList.add('active');
        });
        document.getElementById('btn-compare').disabled = true;
        document.getElementById('input-k').disabled = true;
        document.getElementById('input-m').disabled = true;
        this.updateMesh(null);
        this._render();
    }

    clearMesh() {
        [this.deformedMesh, this.ghostMesh].forEach(m => { if (m) { this.scene.remove(m); m.geometry.dispose(); } });
        this.deformedMesh = this.ghostMesh = null;
    }
}
