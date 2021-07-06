use std::time::Instant;

use futures::executor::block_on;
use tracing::{debug, info, trace};
use tracing_subscriber::EnvFilter;
use winit::{
    event::{Event, WindowEvent},
    event_loop::{ControlFlow, EventLoop},
    window::WindowBuilder,
};

use crate::{
    args::Args,
    graphics::{DrawError, Renderer, Transform},
    input::InputHandler,
    threemf, Mesh, Model,
};

pub fn run_model(model: impl Model) -> anyhow::Result<()> {
    let args = init();

    let mesh = model.instantiate();
    run_inner(mesh, args)?;

    Ok(())
}

pub fn run_mesh(mesh: impl Into<Mesh>) -> anyhow::Result<()> {
    let args = init();
    run_inner(mesh, args)?;
    Ok(())
}

fn init() -> Args {
    tracing_subscriber::fmt()
        .with_env_filter(
            EnvFilter::from_default_env()
                .add_directive("fj=debug".parse().unwrap()),
        )
        .pretty()
        .init();

    Args::parse()
}

fn run_inner(mesh: impl Into<Mesh>, args: Args) -> anyhow::Result<()> {
    info!("Converting geometry to triangle mesh...");

    let start_of_conversion = Instant::now();
    let mesh = mesh.into();
    let conversion_duration = start_of_conversion.elapsed();

    info!(
        "Converted geometry in {}.{:03}s",
        conversion_duration.as_secs(),
        conversion_duration.subsec_millis()
    );

    if let Some(path) = args.export {
        info!("Exporting to `{}`", path.display());
        threemf::export(&mesh, path)?;
        return Ok(());
    }

    debug!("Initializing event loop...");
    let event_loop = EventLoop::new();

    debug!("Initializing window...");
    let window = WindowBuilder::new()
        .with_title("Fornjot")
        .with_maximized(true)
        .with_decorations(true)
        .with_transparent(false)
        .build(&event_loop)
        .unwrap();

    debug!("Initializing event handler...");
    let mut input_handler = InputHandler::new();

    debug!("Initializing transform...");
    let mut transform = Transform::new();

    debug!("Initializing renderer...");
    let mut renderer = block_on(Renderer::new(&window, mesh.into())).unwrap();

    debug!("Finished initialization.");

    event_loop.run(move |event, _, control_flow| {
        trace!("Handling event: {:?}", event);

        match event {
            Event::WindowEvent {
                event: WindowEvent::CloseRequested,
                ..
            } => {
                *control_flow = ControlFlow::Exit;
            }
            Event::WindowEvent {
                event: WindowEvent::Resized(size),
                ..
            } => {
                renderer.handle_resize(size);
            }
            Event::WindowEvent {
                event: WindowEvent::KeyboardInput { input, .. },
                ..
            } => {
                input_handler.handle_keyboard_input(input, control_flow);
            }
            Event::WindowEvent {
                event: WindowEvent::CursorMoved { position, .. },
                ..
            } => {
                input_handler.handle_cursor_moved(position, &mut transform);
            }
            Event::WindowEvent {
                event: WindowEvent::MouseInput { state, button, .. },
                ..
            } => {
                input_handler.handle_mouse_input(state, button);
            }
            Event::WindowEvent {
                event: WindowEvent::MouseWheel { delta, .. },
                ..
            } => {
                input_handler.handle_mouse_wheel(delta, &mut transform);
            }
            Event::MainEventsCleared => {
                window.request_redraw();
            }
            Event::RedrawRequested(_) => {
                match renderer.draw(&transform) {
                    Ok(()) => {}
                    err @ Err(DrawError(wgpu::SwapChainError::Outdated)) => {
                        // I'm getting this from time to time when resizing the
                        // window. It's not catastrophic.
                        println!("Draw error: {:?}", err);
                    }
                    Err(err) => {
                        panic!("Draw error: {:?}", err);
                    }
                }
            }
            _ => {}
        }
    })
}
